% Script to quantify the amount of movement during the eyes open / eyes
% closed task with and without the presence of EEG.
% Zelekha Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all;
addpath('R:\MATLAB\OptiTrack');
% Choose directory to save results in
save_dir = uigetdir('R:\EEG MEG Project\Results','Choose directory to save results in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First load in the OT data
% For all participants:
Pnum = {'P1','P2','P3','P4','P5','P6','P7','P8','P9','P12'}; % OPM-MEG with EEG
% Pnum = {'P2','P3','P4','P5','P6','P7','P8','P9','P10','P12'}; % OPM-MEG alone
for sub = 1:length(Pnum)
% User selection of exported OT file (.csv)
[filename.OptiTrack,path.OptiTrack] = uigetfile('.csv',['Select OT data for ',char(Pnum(sub))],'R:\'); 
opti_data = autoread_opti_v2([path.OptiTrack '\' filename.OptiTrack]); % Load optitrack data
time_init = opti_data.Time;
opti_Fs = length(opti_data.Time)/opti_data.Time(end); % Get sampling frequency

% Plot OT data
helmet_trans_init = [opti_data.Helmet.X_Position,opti_data.Helmet.Y_Position,opti_data.Helmet.Z_Position]*1000; % in mm
% make relative to initial position/orientation
helmet_trans_init = helmet_trans_init-helmet_trans_init(1,:);
figure('color','w'); subplot(1,3,1);
plot(time_init,helmet_trans_init); title('Helmet Translation');
xlabel('Time (s)'); ylabel('Translation (mm)');
xlim([time_init(1) time_init(end)]); legend('x','y','z');
drawnow

% Remove any points where the OT lost the markers
lost_marker = find(abs(helmet_trans_init(:,2)) > 1200);
if lost_marker
    % Chop out data 0.5s either side of the lost markers
    helmet_trans_init(lost_marker(1)-0.5*opti_Fs:lost_marker(end)+0.5*opti_Fs,:) = [];
    time_init(lost_marker(1)-0.5*opti_Fs:lost_marker(end)+0.5*opti_Fs) = [];
end
subplot(1,3,2);
plot(time_init,helmet_trans_init); title('Helmet Translation');
xlabel('Time (s)'); ylabel('Translation (mm)');
xlim([time_init(1) time_init(end)]); legend('x','y','z');
drawnow

skip = input('Type in the number of seconds to skip before the experiment started to exclude non-task-related movements ');
new_init = find(time_init>skip);
helmet_trans = helmet_trans_init;
helmet_trans(1:new_init(1),:) = [];
helmet_trans = helmet_trans - helmet_trans(1,:);
time = time_init;
time(1:new_init(1),:) = [];

subplot(1,3,3); 
plot(time,helmet_trans); title('Helmet Translation');
xlabel('Time (s)'); ylabel('Translation (mm)');
xlim([time(1) time(end)]); legend('x','y','z');
drawnow

pause

max_mvmnt_subs(sub) = max(abs(helmet_trans(:)));
end

cd(save_dir)
save('max_mvmnt_subs.mat','max_mvmnt_subs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get results to see if participants move more with EEG
cd('R:\EEG MEG Project\Results\EEGMEG_MEGdata\Eyes_open_closed');
EEG_MEG_mvmnt = load("max_mvmnt_subs.mat");
EEG_MEG_mvmnt = EEG_MEG_mvmnt.max_mvmnt_subs;
mean(EEG_MEG_mvmnt)
std(EEG_MEG_mvmnt)

cd('R:\EEG MEG Project\Results\MEG_only\Eyes_open_closed');
MEG_mvmnt = load("max_mvmnt_subs.mat");
MEG_mvmnt = MEG_mvmnt.max_mvmnt_subs;
mean(MEG_mvmnt)
std(MEG_mvmnt)


