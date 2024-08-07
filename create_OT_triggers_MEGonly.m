% Script to create markers based on the OptiTrack for the right index
% finger abduction task in MEG. Trigger is when finger returns to
% initial position (i.e. end of finger abduction)
% For this, we need the MEG data and the OT data
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
% should now be called Finger_Label). Export OT data as usual by
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
trig_inds = cat(1,trig_inds,trig_inds1(find((diff(trig_inds1)>max(diff(trig_inds1))/2))+1));

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
[filename.MEG,path.MEG] = uigetfile('.cmeg','Select MEG finger abduction data','M:\'); 
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
save('OT_trigs_MEG.mat','inds')
