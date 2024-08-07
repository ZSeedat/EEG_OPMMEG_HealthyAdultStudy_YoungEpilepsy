% Script to read in EEG data for further analysis using FieldTrip
% Triggers on movement offset using the optitrack markers.
% Bad trials are removed and data are saved out for further analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all; clc;
% Add fieldtrip path
addpath("C:\Users\zseedat\OneDrive - Young Epilepsy\Documents\MATLAB\fieldtrip-20220212\fieldtrip-20220212")
ft_defaults
% Choose directory to save results in
save_dir = uigetdir('R:\EEG MEG Project\Preprocessed_Data\','Choose directory to save processed data in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in the data
% For all participants:
Pnum = {'P1','P2','P3','P5','P6','P7','P8','P9','P10','P11','P12'}; % Exclude bad datasets from this list (e.g. no EEG_MEG sync)
for sub = 1:length(Pnum)
    clearvars -except Pnum sub save_dir
    % User selection of relevant EEG file
    [DataFile,DataPath]=uigetfile('.vhdr',['Select EEG data for ',Pnum{sub}],'R:\');
    % use Fieldtrip
    cfg = [];
    cfg.dataset = [DataPath,DataFile];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49.5 50.5; 99.5 100.5];
    % cfg.dftfilter = 'yes';
    % cfg.dftfreq = [50 100 150];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 300];
    % cfg.reref = 'yes';
    % cfg.refmethod = 'avg';
    % cfg.refchannel = 'all';

    % Read in data
    data_all_chans_trials = ft_preprocessing(cfg);
    EEG_Fs = data_all_chans_trials.hdr.Fs;

    % At the minute all the data is one long trial
    data_all = data_all_chans_trials.trial{1};
    data_all_labels = data_all_chans_trials.label;

    % Remove ECG channel
    ECG_chan = find(strcmp(data_all_chans_trials.label,'ECG'));
    data_all_chans_trials.label(ECG_chan) = [];
    data_all_chans_trials.trial{1}(ECG_chan,:)=[];
    data_all(ECG_chan,:) = [];
    data_all_labels(ECG_chan,:) = [];
    Nchans = size(data_all,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import markers and create trials:
    % % Extract MEG-EEG markers
    % cfg.dataset = [DataPath,DataFile];
    % event = ft_read_event(cfg.dataset,'detectflank',[]);
    % 
    % % Write down the name of the EEG marker you would like to use. Look
    % % for it in the 'event' variable in the MATLAB workspace.
    % event_vals = {event.value};
    % event_vals = char(event_vals{:});
    % sel = ismember(event_vals,{'R128'}); % marker for stim cue
    % trig_inds = find(sel);
    % cue_marker_samp = [];
    % for m = 1:length(trig_inds)
    %     cue_marker_samp(m) = event(trig_inds(m)).sample;
    % end
    % figure; plot(cue_marker_samp,ones(size(cue_marker_samp)),'x')
    % % Clean up triggers:
    % EEG_inds = cue_marker_samp;
    % diff_inds = diff(EEG_inds)>1*EEG_Fs; % trigger indices are more than a second apart
    % diff_inds = cat(2,[0],diff_inds); % Remove EEG-MEG sync at the beginning
    % EEG_inds = EEG_inds(logical(diff_inds));
    % hold on; plot(EEG_inds,ones(size(EEG_inds)),'o')
    % title('Check these are right before proceeding!')
    % drawnow; pause;

    % Load in OptiTrack marker file:
    cd(DataPath)
    save_name = [DataFile(1:end-5),'_OT_velocity_trigs.mat'];
    load([DataPath,save_name]) % trig when finger abduction ends

    % Chop each trial
    offset = -2; % 1 second after cue to make sure they've had time to close their eyes
    seg_duration = 5; % 5s segments to average over
    data = data_all_chans_trials;
    data.trial = [];
    data.time = [];
    for t = 1:length(EEG_inds)
        data.trial{t} = data_all(:,EEG_inds(t)+EEG_Fs*offset:EEG_inds(t)+seg_duration*EEG_Fs+offset*EEG_Fs); % -10ms to 140ms around each marker
        data.time{t} = linspace(offset,seg_duration,size(data.trial{t},2));
        data.sampleinfo(t,:) = [EEG_inds(t)+EEG_Fs*offset, EEG_inds(t)+seg_duration*EEG_Fs+offset*EEG_Fs];
    end
    
    % % Epoch resting state data as it is one long trial at the moment:
    % offset = 0;
    % trial_duration = 7*60; % a single 7 min trial
    % epoch_dur = 5; % each epoch is 5s long
    % epoch = (1:floor(trial_duration/epoch_dur)).*epoch_dur; 
    % data = data_all_chans_trials;
    % if epoch        
    %     data.trial = [];
    %     data.time = [];
    %     for t = 1:length(epoch)
    %         data.trial{t} = data_all(:,EEG_inds(1)+EEG_Fs*offset+(epoch(t)-epoch_dur)*EEG_Fs:EEG_inds(1)+offset*EEG_Fs+epoch(t)*EEG_Fs); % 5s epochs
    %         data.time{t} = linspace(offset,epoch_dur,size(data.trial{t},2));
    %         data.sampleinfo(t,:) = [EEG_inds(1)+EEG_Fs*offset, EEG_inds(1)+epoch_dur*EEG_Fs+offset*EEG_Fs];
    %     end
    % end

    % Check number of trials is correct
    Ntrials = length(data.trial);
    disp([num2str(Ntrials) ' trials found']) % Print to workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visually inspect and reject bad trials
    cfg = [];
    % cfg.method = 'trial';
    cfg.method = 'summary';
    data = ft_rejectvisual(cfg, data);
    good_trials = data.cfg.trials;
    bad_trials = ~(ismember(1:Ntrials,good_trials));
    bad_chans = data.cfg.badchannel;
    % Also remove the bad trials from the trials matrix
    EEG_inds = EEG_inds(good_trials);
    % % Also remove the bad epochs from the trials matrix
    % if epoch
    %     epoch = epoch(good_trials);
    % end
    % % Indices of beginning of good epochs
    % ind1 = EEG_inds(1); 
    % EEG_inds = [];
    % EEG_inds = ind1+EEG_Fs*offset+(epoch-epoch_dur)*EEG_Fs;

    Ntrials = size(good_trials,2);
    disp([num2str(Ntrials) ' trials remain']) % Print to workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Re-reference to the average (without bad channels)
    cfg = [];
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    refchans = data.label; % excludes bad channels
    % %%%% For alpha task only %%%%
    % % Remove any occipital channels
    % refchans(find(contains(refchans,"O"))) = [];
    % refchans(find(startsWith(refchans,"P"))) = [];
    cfg.refchannel = refchans;

    data_reref = ft_preprocessing(cfg, data_all_chans_trials);
    % Over-write data_all to be the re-referenced data
    data_all = data_reref.trial{1};
    data_all_labels = data_reref.label;

    % data_all_labels = data_all_chans_trials.label;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save out the pre-processed data for further analyses
    cd(save_dir)
    mkdir([save_dir,'\',Pnum{sub}])
    cd([save_dir,'\',Pnum{sub}])
    save('data_all.mat','data_all','EEG_Fs','data_all_labels')
    save('good_trial_markers.mat','EEG_inds')
    save('bad_channels.mat','bad_chans')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save out data to view in AnyWave:
% dat = ft_fetch_data(data_reref);
% hdr = ft_fetch_header(data_reref);
% save_path = uigetdir('M:\','Save Folder');
% save_file = inputdlg('Choose a name for your .vhdr file.');
% filename = [save_path,'\',save_file{1},'_EEG.vhdr'];
% ft_write_data(filename, dat, 'header', hdr, 'event', event);
