% Script to create markers based on the OptiTrack for the right index
% finger abduction task in EEG-MEG. Trigger is when finger returns to
% initial position (i.e. end of finger abduction)
% For this, we need the EEG data, MEG data, and OT data
% Zelekha Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all;
addpath('C:\Documents\MATLAB\OptiTrack');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First load in the OT data
% The finger marker should have been labelled in the motive software first.
% Do this by first creating a MarkerSet asset and calling it 'Finger', then
% open the MArkerSet pane and add a marker, then go to the labelling pane,
% select 'Finger' from the drop down list, select QuickLabelling mode, then
% hover over the marker in the Perspective View before clicking on it (it
% should now be called Finger_Label. Export OT data as usual by
% clicking File and selecting Export Tracking Data...

% User selection of exported OT file (.csv)
[filename.OptiTrack,path.OptiTrack] = uigetfile('.csv','Select OT data','M:\'); 
opti_data = autoread_opti_v2([path.OptiTrack '\' filename.OptiTrack]); % Load optitrack data
time_init = opti_data.Time;
opti_Fs = length(opti_data.Time)/opti_data.Time(end); % Get sampling frequency

% Plot OT data
finger_trans_init = [opti_data.Finger_Label.X_Position,opti_data.Finger_Label.Y_Position,opti_data.Finger_Label.Z_Position]*1000; % in mm
% make relative to initial position/orientation
finger_trans_init = finger_trans_init-finger_trans_init(1,:);
figure('color','w');
plot(time_init,finger_trans_init); title('Finger Translation');
xlabel('Time (s)'); ylabel('Translation (mm)');
xlim([time_init(1) time_init(end)]); legend('x','y','z');
drawnow

skip = input('Type in the number of seconds to skip before the experiment started to exclude non-task-related movements ');
new_init = find(time_init>skip);
finger_trans = finger_trans_init;
finger_trans(1:new_init(1),:) = [];
time = time_init;
time(1:new_init(1),:) = [];

% X-translation is most reliable for thresholding
finger_x = abs(finger_trans(:,1));
thresh_val = 1*std(finger_x);
finger_xthresh = double(finger_x > thresh_val);
diff_inds = diff(finger_xthresh);
% When finger returns to initial position
trig_inds1 = find(diff_inds==-1); 

% Clean up trigger indices
trig_inds = trig_inds1(1); % start of first trigger
trig_inds = cat(1,trig_inds,trig_inds1(find((diff(trig_inds1)>max(diff(trig_inds1))/2)+1)));

% Plot check
figure('color','w');
plot(time, finger_xthresh); title('X-translation, thresholded');
xlabel('Time (s)'); ylabel('Finger Translation');
xlim([time(1) time(end)]);
hold on; plot(time(trig_inds), ones(size(trig_inds)).*0.5, 'x')
display(['Number of trials found: ', num2str(length(trig_inds))])

% Correct to start of recording by adding the skip back in
trig_inds_corrected = round(trig_inds + skip*opti_Fs);
figure('color','w');
plot(time_init,finger_trans_init); title('Finger Translation');
xlabel('Time (s)'); ylabel('Translation (mm)');
xlim([time_init(1) time_init(end)]); legend('x','y','z');
hold on; plot(time_init(trig_inds_corrected), ones(size(trig_inds_corrected)).*20, 'x')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put the markers into MEG space:
% User selection of the right index finger abduction MEG file
[filename.MEG,path.MEG] = uigetfile('.cmeg','Select MEG finger abduction data for OT trigger sync','M:\'); 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Put the markers into EEG space:
% User selection of the EEG dataset
[filename.EEG,path.EEG]=uigetfile('.vhdr','Select EEG data','M:\');

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
save_name = [filename.EEG(1:end-5),'_OT_trigs.mat'];
save(save_name,'EEG_inds')

