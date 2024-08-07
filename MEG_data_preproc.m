% Script to read in MEG data for further analysis using FieldTrip.
% Triggers on movement offset using the optitrack markers.
% Bad trials are removed and data are saved out for further analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all; clc;
% Add fieldtrip path
addpath("R:\MATLAB\fieldtrip-20220212\fieldtrip-20220212")
ft_defaults
% Add path with data extraction code in it
addpath('R:\MATLAB')
% Choose directory to save results in
save_dir = uigetdir('R:\EEG MEG Project\Preprocessed_Data\','Choose directory to save processed data in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in the data
% Choose participants:
Pnum = {'P1','P2','P3','P5','P6','P7','P8','P9','P10','P11','P12'}; % Exclude participants where the expeirment failed for whatever reason
for sub = 1:length(Pnum)
    clearvars -except Pnum sub save_dir
    % User selection of the right index finger abduction MEG file
    % Note! Select the .cmeg file with the number at the end (e.g., _001)
    [filename.MEG,path.MEG] = uigetfile('.cmeg',['Select MEG data for ',Pnum{sub}],'R:\');
    MEG_file = [path.MEG,filename.MEG];
    % Put into fieldtrip format
    FT_data_struct = convert_2_fieldtrip(MEG_file); % No preprocessing is done here, the conversion factor is applied and the format is changed

    % Preprocess using Fieldtrip
    cfg = [];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49.5 50.5; 99.5 100.5];
    % cfg.dftfilter = 'yes';
    % cfg.dftfreq = [50 100 150];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 300];

    % Read in data
    data = ft_preprocessing(cfg,FT_data_struct.data_all);
    MEG_Fs = data.fsample;

    % At the minute all the data is one long trial
    data_all = data.trial{1};
    time = data.time{1};
    Nchans = size(data_all,1);

    % Order of sensor positions in helmet:
    sens_order = FT_data_struct.sens_pos_order;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import markers and create trials:
    % Load in the OptiTrack marker file:
    cd(path.MEG)
    load('OT_velocity_trigs_MEG.mat') % trig when finger abduction ends 
    
    % OR

    % % Calculate markers from trigger channel 1
    % trig1 = FT_data_struct.trigger(1,:);
    % figure; plot(trig1)
    % diff_trig1 = diff(trig1);
    % inds = find(diff_trig1 > 3)+1;
    % % Clean up triggers:
    % diff_inds = diff(inds)>1*MEG_Fs; % trigger indices are more than a second apart
    % diff_inds = cat(2,[0],diff_inds); % Remove EEG-MEG sync at the beginning
    % inds = inds(logical(diff_inds));
    % Plot check:
    % hold on; plot(inds,3.*ones(size(inds)),'o');
    % title('Check these are right before proceeding!')
    % drawnow; pause;

    % Chop each trial
    if inds
        offset = -2;
        trial_duration = 5; % a single 7 min trial
        data.trial = [];
        data.time = [];
        for t = 1:length(inds)
            data.trial{t} = data_all(:,inds(t)+MEG_Fs*offset:inds(t)+trial_duration*MEG_Fs+offset*MEG_Fs); % -10ms to 140ms around each marker
            data.time{t} = linspace(offset,trial_duration,size(data.trial{t},2));
            data.sampleinfo(t,:) = [inds(t)+MEG_Fs*offset, inds(t)+trial_duration*MEG_Fs+offset*MEG_Fs];
        end
    end
    
    % % Epoch resting state data as it is one long trial at the moment:
    % offset = 0;
    % trial_duration = 7*60; % a single 7 min trial
    % epoch_dur = 5; % each epoch is 5s long
    % epoch = (1:floor(trial_duration/epoch_dur)).*epoch_dur; 
    % if epoch        
    %     data.trial = [];
    %     data.time = [];
    %     for t = 1:length(epoch)
    %         data.trial{t} = data_all(:,inds(1)+MEG_Fs*offset+(epoch(t)-epoch_dur)*MEG_Fs:inds(1)+offset*MEG_Fs+epoch(t)*MEG_Fs); % 5s epochs
    %         data.time{t} = linspace(offset,epoch_dur,size(data.trial{t},2));
    %         data.sampleinfo(t,:) = [inds(1)+MEG_Fs*offset, inds(1)+epoch_dur*MEG_Fs+offset*MEG_Fs];
    %     end
    % end

    % Check number of trials is correct
    Ntrials = length(data.trial);
    disp([num2str(Ntrials) ' trials found']) % Print to workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visually inspect and reject bad trials / bad channels
    cfg = [];
    cfg.method = 'summary';
    data = ft_rejectvisual(cfg, data);
    good_trials = data.cfg.trials;
    bad_trials = ~(ismember(1:Ntrials,good_trials));
    bad_chans = data.cfg.badchannel;
    % Also remove the bad trials from the trials matrix
    if inds
        inds = inds(good_trials);
    end
    % % Also remove the bad epochs from the trials matrix
    % if epoch
    %     epoch = epoch(good_trials);
    % end
    % % Indices of beginning of good epochs
    % ind1 = inds(1); 
    % inds = [];
    % inds = ind1+MEG_Fs*offset+(epoch-epoch_dur)*MEG_Fs;

    Ntrials = size(good_trials,2);
    disp([num2str(Ntrials) ' trials remain']) % Print to workspace
    data_all_labels = data.cfg.previous.channel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Homogeneous field correction
% As in Tierney et al. 2021, https://www.sciencedirect.com/science/article/pii/S1053811921007576

    % Apply HFC to unfiltered data:
    cfg = [];
    dataHFC = ft_preprocessing(cfg,FT_data_struct.data_all);

    % At the minute all the data is one long trial
    dataHFC_all = dataHFC.trial{1};

    % Remove bad channels before mean field correcting
    data_labels_HFC = string(data_all_labels);
    bad_chan_inds = [];
    if length(bad_chans)>0
        for ch = 1:length(bad_chans)
            bad_chan_inds(ch) = find(string(data_all_labels)==bad_chans{ch});
        end
        dataHFC_all(bad_chan_inds,:) = [];
        data_labels_HFC(bad_chan_inds)=[];
    end
    dataHFC_all = dataHFC_all-repmat(mean(dataHFC_all,2),1,size(dataHFC_all,2));
    hfc_check = figure; subplot(2,1,1); 
    for ch = 1:size(dataHFC_all,1)
        p = plot(time,dataHFC_all(ch,:));
        p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow(data_labels_HFC(ch),[]);
        hold on
    end
    title('Data without HFC'); xlabel('Time,s'); ylabel('Field,nT')

    % Load helmet info
    Helm_config = tdfread([MEG_file(1:54),'HelmConfig.tsv'],'\t');
    helm_chan_order = cellstr(Helm_config.Sensor);
    helm_chan_order = helm_chan_order(1:end-1); % the last entry isn't a sensor
    % Add the square brackets in to match the channel labels from '_channels.tsv'
    for n = 1:length(helm_chan_order)
        helm_chan_order{n} = [helm_chan_order{n}(1:3),'[',helm_chan_order{n}(4),']'];
    end

    % First find orientations of each sensor
    sens_pos = ones(length(data_labels_HFC),3);
    sens_ors = ones(length(data_labels_HFC),3);
    for n = 1:length(data_labels_HFC)
        sens_helmet_ind = find(data_labels_HFC(n)==helm_chan_order);
        sens_pos(n,1) = str2double(Helm_config.Px(sens_helmet_ind,:));
        sens_pos(n,2) = str2double(Helm_config.Py(sens_helmet_ind,:));
        sens_pos(n,3) = Helm_config.Pz(sens_helmet_ind,:); % z is already a double
        sens_ors(n,1) = Helm_config.Ox(sens_helmet_ind,:);
        sens_ors(n,2) = Helm_config.Oy(sens_helmet_ind,:);
        sens_ors(n,3) = str2double(Helm_config.Oz(sens_helmet_ind,:)); % string array
    end

    % Check the orientations are right
    figure
    quiver3(sens_pos(:,1),sens_pos(:,2),sens_pos(:,3),...
        sens_ors(:,1),sens_ors(:,2),sens_ors(:,3))
    axis equal

    % Apply homogeneous field correction
    disp('Applying homogeneous field correction');
    N = sens_ors;
    M = eye(length(N)) - N*pinv(N);
    Sdata = M*dataHFC_all;
    figure(hfc_check); subplot(2,1,2); 
    for ch = 1:size(Sdata,1)
        p = plot(time,Sdata(ch,:));
        p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow(data_labels_HFC(ch),[]);
        hold on
    end
    title('Data with HFC'); xlabel('Time,s'); ylabel('Field,nT')

    % Put back into a fieldtrip data format:
    dataHFC.trial{1} = Sdata;
    dataHFC.label = cellstr(data_labels_HFC);
                
    % Filter the HFC data using FieldTrip
    cfg = [];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49.5 50.5; 99.5 100.5];
    % cfg.dftfilter = 'yes';
    % cfg.dftfreq = [50 100 150];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 300];
    dataHFC_f = ft_preprocessing(cfg,dataHFC);
    
    % Convert into 1 long trial to save out
    clear data_all data_all_labels
    data_all = dataHFC_f.trial{1};
    data_all_labels = dataHFC_f.label;

    drawnow; pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save out the pre-processed data for further analyses
    cd(save_dir)
    mkdir([save_dir,'\',Pnum{sub}])
    cd([save_dir,'\',Pnum{sub}])
    save('data_all.mat','data_all','MEG_Fs','data_all_labels')
    save('good_trial_markers.mat','inds')
    save('bad_channels.mat','bad_chans','bad_chan_inds')
    save('sens_order.mat','sens_order')
end


