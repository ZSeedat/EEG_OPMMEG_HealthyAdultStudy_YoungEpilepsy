% Script to read in and plot exported OptiTrack data 
% Also adds in triggers based on velocity to look for movement related
% artefact.
% Before running this script, please make sure all your markers are
% labelled and all your rigid bodies are resolved in Motive. Then export 
% trackiong data as a .csv file which can be read into this script.
% Script from Molly Rea at UoN and edited by Zelekha Seedat at YE
% Last edited: 16/01/2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User selection of OptiTrack dataset exported from Motive
clear; clc;
addpath('R:\MATLAB\OptiTrack')
addpath('R:\MATLAB\fieldtrip-20220212\fieldtrip-20220212')

Pnum = {'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12'};
for sub = 1:length(Pnum)
[filename.OptiTrack,path.OptiTrack] = uigetfile('*.csv',['Select OT finger abduction data for ',Pnum{sub}],'R:\');
opti_data = autoread_opti_v2([path.OptiTrack '\' filename.OptiTrack]); % Load optitrack data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract and plot helmet data
% Organise rotation/translation data
rots1 = [opti_data.Helmet.X_Rotation,opti_data.Helmet.Y_Rotation,opti_data.Helmet.Z_Rotation]; % in degrees
trans1 = [opti_data.Helmet.X_Position,opti_data.Helmet.Y_Position,opti_data.Helmet.Z_Position]*1000; % in mm
time_init = opti_data.Time;
opti_Fs = length(opti_data.Time)/opti_data.Time(end); % Get sampling frequency

% Exclude large head movement at beginnning of recording
% skip = input('Type in the number of seconds to skip before the experiment started to exclude non-task-related movements ');
skip = 0;
new_init = find(time_init>skip);
rots = rots1;
rots(1:new_init(1),:) = [];
trans = trans1;
trans(1:new_init(1),:) = [];
time = time_init;
time(1:new_init(1),:) = [];
% Make relative to initial position/orientation
rots = rots-rots(1,:);
trans = trans-trans(1,:);

% Plot data
figure('color','w'); 
% Rigid body rotations
subplot(2,1,1);
plot(time,rots); title('Helmet Rotation');
xlabel('Time (s)'); ylabel('Rotation (\circ)');
xlim([time(1) time(end)]);
% Rigid body translations
subplot(2,1,2);
plot(time,trans); title('Helmet Translation');
xlabel('Time (s)'); ylabel('Translation (mm)');
xlim([time(1) time(end)]);

% Get max rotation/translation and print to command window
max_r = max(abs(rots)); % x y z
max_t = max(abs(trans)); % x y z
disp(['Maximum helmet rotation about x: ',num2str(max_r(1),3),' degrees']);
disp(['Maximum helmet rotation about y: ',num2str(max_r(2),3),' degrees']);
disp(['Maximum helmet rotation about z: ',num2str(max_r(3),3),' degrees']);
disp(['Maximum helmet translation about x: ',num2str(max_t(1),3),' mm']);
disp(['Maximum helmet translation about y: ',num2str(max_t(2),3),' mm']);
disp(['Maximum helmet translation about z: ',num2str(max_t(3),3),' mm']);

helmet_velocity = diff(trans);
vel_mag = sqrt(helmet_velocity(:,1).^2+helmet_velocity(:,2).^2+helmet_velocity(:,3).^2);
figure('color','w');
plot(time(1:end-1),helmet_velocity)
xlabel('Time,s'); ylabel('Speed, m/s')
hold on; plot(time(1:end-1),vel_mag)

thresh_val = 0.1; %m/s
vel_thresh = double(vel_mag > thresh_val);
diff_inds = diff(vel_thresh);
% When finger returns to initial position
diff_inds = find(diff_inds==1); 

% Clean up trigger indices
trig_inds = diff_inds(1); % start of first trigger
trig_inds = cat(1,trig_inds,diff_inds(find(diff(diff_inds)>3*opti_Fs)+1));

% Plot check
hold on; plot(time(trig_inds), ones(size(trig_inds)).*0.5, 'x')
display(['Number of trials found: ', num2str(length(trig_inds))])

% Correct to start of recording by adding the skip back in
trig_inds_corrected = round(trig_inds + skip*opti_Fs);
helmet_velocity1 = diff(trans1);
vel_mag_init = sqrt(helmet_velocity1(:,1).^2+helmet_velocity1(:,2).^2+helmet_velocity1(:,3).^2);
figure('color','w');
plot(time_init(1:end-1),vel_mag_init,'k-'); title('Helmet velocity');
xlabel('Time (s)'); ylabel({'Velocity','magntiude (m/s)'});
xlim([time_init(1) time_init(end-1)]); 
hold on; plot(time_init(trig_inds_corrected), ones(size(trig_inds_corrected)).*0.1, 'x','MarkerSize',20,'LineWidth',3)
set(gca,'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put the markers into MEG space:
% User selection of the right index finger abduction MEG file
[filename.MEG,path.MEG] = uigetfile('.cmeg','Select MEG finger abduction data','R:\'); 
MEG_data = read_cMEG_data_split([path.MEG,filename.MEG]);
MEG_Fs = MEG_data.samp_frequency;

% Plot trigger channels so you can select the one with the OptiTrack sync
chan2lookat = 1:8; % First 8 channels are triggers
MEG_trigs = MEG_data.data(chan2lookat,:);
figure('Color','w')
plot(MEG_data.time,MEG_trigs)
xlabel('Time,s'); ylabel('Voltage'); title('Trigger Channels')
legend({'chan 1','chan 2','chan 3','chan 4','chan 5','chan 6','chan 7','chan 8'})
set(gca,'FontSize',16)

% Get the optitrack start marker:
OT_start = find(diff(MEG_trigs(2,:)>1)==-1);
OT_start = OT_start(1);

% Now convert the OT indices to match the MEG sample frequency:
MEG_inds = (trig_inds_corrected./opti_Fs).*MEG_Fs;

% And add in the delay between the MEG rec starting and OT starting
MEG_inds = MEG_inds+OT_start;

% Plot check (should be after the trigger cue for each trial)
hold on; plot(MEG_inds/MEG_Fs, ones(size(MEG_inds)),'*')

% Save out for use in subsequent scripts and in analysis app
cd(path.MEG)
inds = MEG_inds;
save('OT_velocity_trigs_MEG.mat','inds')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put the markers into EEG space:
% User selection of the EEG dataset
[filename.EEG,path.EEG]=uigetfile('.vhdr','Select EEG data','R:\');

% Use Fieldtrip to read in data
cfg = [];
cfg.dataset = [path.EEG,filename.EEG];
EEG_data = ft_preprocessing(cfg);
EEG_Fs = EEG_data.hdr.Fs;

% Convert the OT indices to match the EEG sample frequency:
EEG_inds = (trig_inds_corrected./opti_Fs).*EEG_Fs;

% Find the EEG/MEG sync in the EEG data
EEG_event = ft_read_event(cfg.dataset,'detectflank',[]);
event_vals = {EEG_event.value};
event_vals = char(event_vals{:});
sel = ismember(event_vals,{'R128'});
EEG_trig_inds = find(sel);
EEG_MEG_sync = EEG_event(EEG_trig_inds(1)).sample; % in the EEG data

% Find the MEG/EEG sync in the MEG data
trig1_high = find(diff(MEG_trigs(1,:)>1)==1); 
MEG_EEG_sync = trig1_high(1); % in the MEG data

% Find the delay between the MEG_EEG trigger and the OT start
OT_MEGsync_delay = (OT_start-MEG_EEG_sync)/MEG_Fs; % In seconds
OT_EEGsync_delay = ((OT_start-MEG_EEG_sync)/MEG_Fs)*EEG_Fs; % In EEG samples

% And add in the delay between the EEG rec starting and OT starting
EEG_inds = round(EEG_inds)+EEG_MEG_sync+round(OT_EEGsync_delay);

% Plot check (should be after the trigger cue for each trial)
% Get the EEG cue markers
for i = 1:length(EEG_trig_inds) 
    EEG_cue_inds(i) = EEG_event(EEG_trig_inds(i)).sample;
end
% discount the 1st one which is the trig sync
EEG_cue_inds(1) = [];
EEG_pseudo_trig_chan = zeros(size(EEG_data.time{1}));
EEG_pseudo_trig_chan(EEG_cue_inds) = 1;
EEG_pseudo_trig_chan(EEG_MEG_sync) = 1;
figure; plot(EEG_data.time{1},EEG_pseudo_trig_chan); hold on;
plot(EEG_data.time{1}(EEG_inds),ones(length(EEG_inds),1).*0.5,'*')

% Save out for use in subsequent scripts and in analysis app
cd(path.EEG)
save_name = [filename.EEG(1:end-5),'_OT_velocity_trigs.mat'];
save(save_name,'EEG_inds')

end

