% Script to read in MEG data for further analysis using FieldTrip
% Gets TFS and beta amplitude timecourse for highest SNR channel per
% subject and then averages over subjects
% Triggers on movement offset using the optitrack markers
% Zelekha Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all; clc;
% Add fieldtrip path
addpath("C:\Documents\MATLAB\fieldtrip-20220212\fieldtrip-20220212")
ft_defaults
% Add path with data extraction code in it
addpath('C:\Documents\MATLAB')
% Choose directory to save results in
save_dir = uigetdir('M:\EEG MEG Project\Results','Choose directory to save results in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in the data
% For all participants:
Pnum = {'P1','P2','P3','P5','P6','P7','P10','P11','P12'};
for sub = 1:length(Pnum)
    clearvars -except Pnum sub save_dir
    % User selection of the right index finger abduction MEG file
    [filename.MEG,path.MEG] = uigetfile('.cmeg',['Select MEG finger abduction data for ',Pnum{sub}],'M:\');
    MEG_file = [path.MEG,filename.MEG];
    % Put into fieldtrip format
    FT_data_struct = convert_2_fieldtrip(MEG_file);

    % Preprocess using Fieldtrip
    cfg = [];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49.5 50.5];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 150];

    % Read in data
    data = ft_preprocessing(cfg,FT_data_struct.data_all);
    MEG_Fs = data.fsample;

    % At the minute all the data is one long trial
    data_all = data.trial{1};
    Nchans = size(data_all,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import markers and create trials:
    % Load in the OptiTrack marker file:
    cd(path.MEG)
    load('OT_trigs_MEG.mat') % trig when finger abduction ends

    % Chop each trial
    offset = -1;
    seg_duration = 5; % 5s segments to average over
    data.trial = [];
    data.time = [];
    for t = 1:length(inds)
        data.trial{t} = data_all(:,inds(t)+MEG_Fs*offset:inds(t)+seg_duration*MEG_Fs+offset*MEG_Fs); % -10ms to 140ms around each marker
        data.time{t} = linspace(offset,seg_duration,size(data.trial{t},2));
        data.sampleinfo(t,:) = [inds(t)+MEG_Fs*offset, inds(t)+seg_duration*MEG_Fs+offset*MEG_Fs];
    end

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
    % Also remove the bad trials from the trials matrix
    inds = inds(good_trials);
    Ntrials = size(good_trials,2);
    disp([num2str(Ntrials) ' trials remain']) % Print to workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot time-frequency spectrogram
    highpass = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 65 70];% 75 80 85 90 95 100 105 110];
    lowpass = [4 6 8 10 13 20 25 30 35 40 45 50 55 60 65 70 75 80];% 85 90 95 100 105 110 115 120];
    fre = highpass + ((lowpass - highpass)./2);
    % Control window
    % duration = input('Please type in the trial duration in seconds, then press enter.');
    duration = 5; % 5 second trials
    MRBD_win = ([-1, 0]-offset).*MEG_Fs+1; % desync window is 1 to 0s before movement cessation
    con_win = [duration-0.5 duration].*MEG_Fs;
    trl_time_induced = linspace(offset,duration+offset,duration*MEG_Fs);

    % Load in sensor layout for plotting
    load('M:\Helmet_info\purple_adult_L.mat') % the same helmet was used for everyone
    layout = Helmet_info;
    lay = layout.lay;

    % Initialise figures for plotting
    h = figure(200);
    if exist('lay')
        p = uipanel('Title','Axes','Units','normalized','Position',...
            [0.85 0.85 0.1 0.12],'BackgroundColor','w','BorderType','none');
        ax = axes(p);
        xlabel('Time (s)');ylabel('Frequency (Hz)')
        ylim([min(fre) max(fre)])
        xlim([min(trl_time_induced) max(trl_time_induced)])
        axis on
        colorbar;caxis([-0.5 0.5])
    end
    h = figure(300);
    if exist('lay')
        p = uipanel('Title','Axes','Units','normalized','Position',...
            [0.85 0.85 0.1 0.12],'BackgroundColor','w','BorderType','none');
        ax = axes(p);
        xlabel('Time (s)');ylabel('Frequency (Hz)')
        ylim([min(fre) max(fre)])
        xlim([min(trl_time_induced) max(trl_time_induced)])
        axis on
        colorbar;caxis([-0.5 0.5])
    end

    % Get the indices for Y and Z channels
    Y_inds = 1:64;
    Z_inds = 65:128;

    % Loop over each sensor
    TFS_Y = [];
    TFS_Z = [];
    data_all_labels = data.cfg.previous.channel;
    Nchans = length(Y_inds);
    for n = 1:Nchans
        % Y-channels first
        MEG_ch_Y = data_all(Y_inds(n),:)';
        % Filter data within bands and calculate envelope
        MEG_ch_fb_Y = zeros(length(MEG_ch_Y),length(fre));
        for fb = 1:length(highpass)
            fprintf('\n Band %u/%u ',fb,length(highpass))
            filt_dat_Y = nut_filter3_nottsapp(MEG_ch_Y,'butter','bp',3,highpass(fb),lowpass(fb),MEG_Fs,1);
            MEG_ch_fb_Y(:,fb) = abs(hilbert(filt_dat_Y));
        end
        MEG_mean_Y = zeros(duration*MEG_Fs,length(fre));
        % Chop data into trials
        for fb = 1:length(highpass)
            MEG_ch_fb_trials_Y = [];
            for i = 1:Ntrials
                MEG_ch_fb_trials_Y = cat(1,MEG_ch_fb_trials_Y,MEG_ch_fb_Y(inds(i)+offset*MEG_Fs+1:inds(i)+offset*MEG_Fs+(duration*MEG_Fs),fb));
            end
            MEG_ch_filt_Y = reshape(MEG_ch_fb_trials_Y,duration*MEG_Fs,Ntrials);
            % Average across trials
            MEG_mean_Y(:,fb) = mean(MEG_ch_filt_Y,2);
        end
        meanrest_Y = mean(MEG_mean_Y(con_win(1):con_win(2),:),1);
        meanrestmat_Y = repmat(meanrest_Y,size(MEG_mean_Y,1),1);
        TFS_Y = (MEG_mean_Y'-meanrestmat_Y')./meanrestmat_Y';

        % Z-channels
        MEG_ch_Z = data_all(Z_inds(n),:)';
        % Filter data within bands and calculate envelope
        MEG_ch_fb_Z = zeros(length(MEG_ch_Z),length(fre));
        for fb = 1:length(highpass)
            fprintf('\n Band %u/%u ',fb,length(highpass))
            filt_dat_Z = nut_filter3_nottsapp(MEG_ch_Z,'butter','bp',3,highpass(fb),lowpass(fb),MEG_Fs,1);
            MEG_ch_fb_Z(:,fb) = abs(hilbert(filt_dat_Z));
        end
        MEG_mean_Z = zeros(duration*MEG_Fs,length(fre));
        % Chop data into trials
        for fb = 1:length(highpass)
            MEG_ch_fb_trials_Z = [];
            for i = 1:Ntrials
                MEG_ch_fb_trials_Z = cat(1,MEG_ch_fb_trials_Z,MEG_ch_fb_Z(inds(i)+offset*MEG_Fs+1:inds(i)+offset*MEG_Fs+(duration*MEG_Fs),fb));
            end
            MEG_ch_filt_Z = reshape(MEG_ch_fb_trials_Z,duration*MEG_Fs,Ntrials);
            % Average across trials
            MEG_mean_Z(:,fb) = mean(MEG_ch_filt_Z,2);
        end
        meanrest_Z = mean(MEG_mean_Z(con_win(1):con_win(2),:),1);
        meanrestmat_Z = repmat(meanrest_Z,size(MEG_mean_Z,1),1);
        TFS_Z = (MEG_mean_Z'-meanrestmat_Z')./meanrestmat_Z';

        % Plot
        % Y-channel plot
        figure(200)
        if exist('lay')
            p = uipanel('Title',lay.label(n),...
                'Units','normalized','Position',...
                [(0.5+lay.pos(n,1))*0.95 (0.5+lay.pos(n,2))*0.95 0.07 0.09],...
                'BackgroundColor','w','BorderType','none');
            ax = axes(p);
        else
            subplot(ceil(sqrt(Nchans)),ceil(sqrt(Nchans)),n);
        end
        pcolor(trl_time_induced,fre,TFS_Y);shading interp
        axis fill
        caxis([-0.5 0.5])
        set(gcf,'color',[1 1 1])
        if ~exist('lay')
            xlabel('Time (s)');ylabel('Frequency (Hz)')
            title(lay.label(n))
        end
        axis off
        title(sprintf(data_all_labels{Y_inds(n)}))
        pause(0.1)
        
        % Z-channel plot
        figure(300)
        if exist('lay')
            p = uipanel('Title',lay.label(n),...
                'Units','normalized','Position',...
                [(0.5+lay.pos(n,1))*0.95 (0.5+lay.pos(n,2))*0.95 0.07 0.09],...
                'BackgroundColor','w','BorderType','none');
            ax = axes(p);

        else
            subplot(ceil(sqrt(64)),ceil(sqrt(64)),n);
        end
        pcolor(trl_time_induced,fre,TFS_Z);shading interp
        axis fill
        caxis([-0.5 0.5])
        set(gcf,'color',[1 1 1])
        if ~exist('lay')
            xlabel('Time (s)');ylabel('Frequency (Hz)')
            title(lay.label(n))
        end
        axis off
        title(sprintf(data_all_labels{Z_inds(n)}))
        pause(0.1)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select the sensor with the greatest beta-band SNR to save out
    beta_highpass = 13;
    beta_lowpass = 30;
    PMBR_win = ([1, 2]-offset).*MEG_Fs; % rebound window is 1 to 2s after movement cessation
    MRBD_win = ([-1, 0]-offset).*MEG_Fs+1; % desync window is 1 to 0s before movement cessation
    con_win = [duration-0.5 duration].*MEG_Fs;
    for n = 1:Nchans
        % Y-channels
        MEG_ch_Y = data_all(Y_inds(n),:)';
        % Filter data into beta bands and calculate envelope
        filt_dat_Y = nut_filter3_nottsapp(MEG_ch_Y,'butter','bp',3,beta_highpass,beta_lowpass,MEG_Fs,1);
        MEG_ch_beta_Y = abs(hilbert(filt_dat_Y));
        % Chop data into trials
        MEG_ch_beta_trials_Y = [];
        for i = 1:Ntrials
            MEG_ch_beta_trials_Y = cat(1,MEG_ch_beta_trials_Y,MEG_ch_beta_Y(inds(i)+offset*MEG_Fs+1:inds(i)+offset*MEG_Fs+(duration*MEG_Fs)));
        end
        MEG_ch_beta_filt_Y = reshape(MEG_ch_beta_trials_Y,duration*MEG_Fs,Ntrials);
        % Average across trials
        MEG_mean_Y = mean(MEG_ch_beta_filt_Y,2);
        meanrest_Y = mean(MEG_mean_Y(con_win(1):con_win(2)));
        meanrestmat_Y = repmat(meanrest_Y,size(MEG_mean_Y,1),1);
        beta_trace_Y(n,:) = (MEG_mean_Y'-meanrestmat_Y')./meanrestmat_Y';
        beta_chan_SNR_Y(n) = (mean(beta_trace_Y(n,PMBR_win(1):PMBR_win(2)))-mean(beta_trace_Y(n,MRBD_win(1):MRBD_win(2))))/std(beta_trace_Y(n,MRBD_win(1):MRBD_win(2)));

        % Z-channels
        MEG_ch_Z = data_all(Z_inds(n),:)';
        % Filter data into beta bands and calculate envelope
        filt_dat_Z = nut_filter3_nottsapp(MEG_ch_Z,'butter','bp',3,beta_highpass,beta_lowpass,MEG_Fs,1);
        MEG_ch_beta_Z = abs(hilbert(filt_dat_Z));
        % Chop data into trials
        MEG_ch_beta_trials_Z = [];
        for i = 1:Ntrials
            MEG_ch_beta_trials_Z = cat(1,MEG_ch_beta_trials_Z,MEG_ch_beta_Z(inds(i)+offset*MEG_Fs+1:inds(i)+offset*MEG_Fs+(duration*MEG_Fs)));
        end
        MEG_ch_beta_filt_Z = reshape(MEG_ch_beta_trials_Z,duration*MEG_Fs,Ntrials);
        % Average across trials
        MEG_mean_Z = mean(MEG_ch_beta_filt_Z,2);
        meanrest_Z = mean(MEG_mean_Z(con_win(1):con_win(2)));
        meanrestmat_Z = repmat(meanrest_Z,size(MEG_mean_Z,1),1);
        beta_trace_Z(n,:) = (MEG_mean_Z'-meanrestmat_Z')./meanrestmat_Z';
        beta_chan_SNR_Z(n) = (mean(beta_trace_Z(n,PMBR_win(1):PMBR_win(2)))-mean(beta_trace_Z(n,MRBD_win(1):MRBD_win(2))))/std(beta_trace_Z(n,MRBD_win(1):MRBD_win(2)));

        % plot
%         figure(400)
%         if exist('lay')
%             p = uipanel('Title',lay.label(n),...
%                 'Units','normalized','Position',...
%                 [(0.5+lay.pos(n,1))*0.95 (0.5+lay.pos(n,2))*0.95 0.07 0.09],...
%                 'BackgroundColor','w','BorderType','none');
%             ax = axes(p);
% 
%         else
%             subplot(ceil(sqrt(Nchans)),ceil(sqrt(Nchans)),n);
%         end
%         plot(trl_time_induced,beta_trace_Y(n,:));
%         axis fill
%         set(gcf,'color',[1 1 1])
%         if ~exist('lay')
%             xlabel('Time (s)');ylabel('Frequency (Hz)')
%             title(lay.label(n))
%         end
%         axis off
%         title(data_all_labels{Y_inds(n)})
%         ylim([-1 1])
%         pause(0.1)

%         figure(500)
%         if exist('lay')
%             p = uipanel('Title',lay.label(n),...
%                 'Units','normalized','Position',...
%                 [(0.5+lay.pos(n,1))*0.95 (0.5+lay.pos(n,2))*0.95 0.07 0.09],...
%                 'BackgroundColor','w','BorderType','none');
%             ax = axes(p);
% 
%         else
%             subplot(ceil(sqrt(Nchans)),ceil(sqrt(Nchans)),n);
%         end
%         plot(trl_time_induced,beta_trace_Z(n,:));
%         axis fill
%         set(gcf,'color',[1 1 1])
%         if ~exist('lay')
%             xlabel('Time (s)');ylabel('Frequency (Hz)')
%             title(lay.label(n))
%         end
%         axis off
%         title(data_all_labels{Z_inds(n)})
%         ylim([-1 1])
%         pause(0.1)
    end
    % Sensor with max SNR in the left centro-parietal electrodes
    relevant_sensors = [8:16,72:80]; % 8th sensor to the 16th sensor, both y and z-channels
    beta_chan_SNR_all = [beta_chan_SNR_Y, beta_chan_SNR_Z];
    [max_SNR,max_SNR_idx] = max(beta_chan_SNR_all(relevant_sensors));
    chan = data_all_labels{relevant_sensors(max_SNR_idx)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot TFS and beta timecourse for selected sensor
    n = find(strcmp(chan,data_all_labels));
    MRBD_win = ([-1, 0]-offset).*MEG_Fs+1; % desync window is 1 to 0s before movement cessation
    con_win = [duration-0.5 duration].*MEG_Fs;
    offset = -1;
    MEG_ch = data_all(n,:)';
    % Filter data within bands and calculate envelope
    MEG_ch_fb = zeros(length(MEG_ch),length(fre));
    for fb = 1:length(highpass)
        fprintf('\n Band %u/%u ',fb,length(highpass))
        filt_dat = nut_filter3_nottsapp(MEG_ch,'butter','bp',3,highpass(fb),lowpass(fb),MEG_Fs,1);
        MEG_ch_fb(:,fb) = abs(hilbert(filt_dat));
    end
    MEG_mean = zeros(duration*MEG_Fs,length(fre));
    % Chop data into trials
    for fb = 1:length(highpass)
        MEG_ch_fb_trials = [];
        for i = 1:Ntrials
            MEG_ch_fb_trials = cat(1,MEG_ch_fb_trials,MEG_ch_fb(inds(i)+offset*MEG_Fs+1:inds(i)+offset*MEG_Fs+(duration*MEG_Fs),fb));
        end
        MEG_ch_filt = reshape(MEG_ch_fb_trials,duration*MEG_Fs,Ntrials);
        % Average across trials
        MEG_mean(:,fb) = mean(MEG_ch_filt,2);
    end
    meanrest = mean(MEG_mean(con_win(1):con_win(2),:),1);
    meanrestmat = repmat(meanrest,size(MEG_mean,1),1);
    TFS = (MEG_mean'-meanrestmat')./meanrestmat';
    figure; subplot(1,2,1)
    pcolor(trl_time_induced,fre,TFS);shading interp
    set(gcf,'color',[1 1 1])
    title(chan)
    pause(0.1)
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    ylim([min(fre) max(fre)])
    xlim([min(trl_time_induced) max(trl_time_induced)])
    axis on
    colorbar; caxis([-1 1])
    set(gca,'FontSize',16)

    % Filter data into beta bands and calculate envelope
    filt_dat = nut_filter3_nottsapp(MEG_ch,'butter','bp',3,beta_highpass,beta_lowpass,MEG_Fs,1);
    MEG_ch_beta = abs(hilbert(filt_dat));

    % Chop data into trials
    MEG_ch_beta_trials = [];
    for i = 1:Ntrials
        MEG_ch_beta_trials = cat(1,MEG_ch_beta_trials,MEG_ch_beta(inds(i)+offset*MEG_Fs+1:inds(i)+offset*MEG_Fs+(duration*MEG_Fs)));
    end
    MEG_ch_beta_filt = reshape(MEG_ch_beta_trials,duration*MEG_Fs,Ntrials);
    % Average across trials
    MEG_mean = mean(MEG_ch_beta_filt,2);

    meanrest = mean(MEG_mean(con_win(1):con_win(2)));
    meanrestmat = repmat(meanrest,size(MEG_mean,1),1);
    beta_trace = (MEG_mean'-meanrestmat')./meanrestmat';
    beta_chan_SNR = (mean(beta_trace(PMBR_win(1):PMBR_win(2)))-mean(beta_trace(MRBD_win(1):MRBD_win(2))))/std(beta_trace(MRBD_win(1):MRBD_win(2)));
    p2p = max(beta_trace(:))-min(beta_trace(:));
    subplot(1,2,2);
    plot(trl_time_induced,beta_trace(:));
    set(gcf,'color',[1 1 1])
    title([chan,' SNR ',num2str(beta_chan_SNR)])
    set(gca,'FontSize',16)
    xlabel('Time (s)'); ylabel('Beta power change from baseline')
    xlim([min(trl_time_induced) max(trl_time_induced)])

    % Save out participant TFS, beta timecourse and SNR
%     cd(save_dir)
%     mkdir([save_dir,'\',Pnum{sub}])
%     cd([save_dir,'\',Pnum{sub}])
%     save('TFS.mat','TFS','trl_time_induced','fre','chan','Ntrials')
%     save('beta_trace.mat','beta_trace','trl_time_induced','chan','Ntrials')
%     save('beta_chan_SNR.mat','beta_chan_SNR','chan','Ntrials')

    drawnow
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot group average result:
% EEG-MEG
% For all participants:
Pnum = {'P1','P2','P3','P5','P6','P7','P10','P11','P12'}; % P4 has no EEG_MEG sync, P8 and P9 had no rebound
save_dir = 'R:\EEG MEG Project\Results\EEGMEG_MEGdata\R_index_abduction';
% clearvars -except Pnum save_dir
for sub = 1:length(Pnum)
    cd([save_dir,'\',Pnum{sub}])
    load("TFS.mat")
    TFS_sub(:,:,sub) = TFS;
    load("beta_trace.mat")
    beta_trace_sub(sub,:) = beta_trace;
    load("beta_chan_SNR.mat")
    beta_SNR_sub_EEGMEG(sub) = beta_chan_SNR;
    Ntrials_sub_EEGMEG(sub) = Ntrials;
end
% Plot average TFS
TFS_avg = mean(TFS_sub,3);
figure; subplot(1,2,1)
pcolor(trl_time_induced,fre,TFS_avg);shading interp
set(gcf,'color',[1 1 1])
pause(0.1)
xlabel('Time (s)');ylabel('Frequency (Hz)')
ylim([min(fre) max(fre)])
xlim([min(trl_time_induced) max(trl_time_induced)])
axis on
colorbar; caxis([-0.5 0.5])
set(gca,'FontSize',16)
% Plot average beta timecourse
beta_trace_avg = mean(beta_trace_sub,1);
beta_trace_std = std(beta_trace_sub,1);
lo = beta_trace_avg - beta_trace_std;
hi = beta_trace_avg + beta_trace_std;
% Define patch
subplot(1,2,2);
hp1 = patch([trl_time_induced, trl_time_induced(end:-1:1)], [lo, hi(end:-1:1)], 'r'); hold on;
h1 = plot(trl_time_induced,beta_trace_avg);
set(hp1, 'facecolor', [173 216 230]/255, 'edgecolor', 'none');
set(h1, 'color', [65 105 225]/255);
set(gcf,'color',[1 1 1])
set(gca,'FontSize',16)
xlabel('Time (s)'); ylabel('Beta power change from baseline')
xlim([min(trl_time_induced) max(trl_time_induced)])
mean_SNR = mean(beta_SNR_sub_EEGMEG);
std_SNR = std(beta_SNR_sub_EEGMEG);
mean_SNR_EEGMEG = mean_SNR;
std_SNR_EEGMEG = std_SNR;
beta_trace_avg_EEGMEG = beta_trace_avg;
beta_trace_std_EEGMEG = beta_trace_std;

% MEG only
% For all participants:
Pnum = {'P1','P2','P3','P4','P5','P6','P7','P10','P12'}; % P11 had no OptiTrack data, P8 and P9 had no rebound
save_dir = 'R:\EEG MEG Project\Results\MEG_only\R_index_abduction';
% clearvars -except Pnum save_dir
for sub = 1:length(Pnum)
    cd([save_dir,'\',Pnum{sub}])
    load("TFS.mat")
    TFS_sub(:,:,sub) = TFS;
    load("beta_trace.mat")
    beta_trace_sub(sub,:) = beta_trace;
    load("beta_chan_SNR.mat")
    beta_SNR_sub_MEGonly(sub) = beta_chan_SNR;
    Ntrials_sub_MEGonly(sub) = Ntrials;
end
% Plot average TFS
TFS_avg = mean(TFS_sub,3);
figure; subplot(1,2,1)
pcolor(trl_time_induced,fre,TFS_avg);shading interp
set(gcf,'color',[1 1 1])
pause(0.1)
xlabel('Time (s)');ylabel('Frequency (Hz)')
ylim([min(fre) max(fre)])
xlim([min(trl_time_induced) max(trl_time_induced)])
axis on
colorbar; caxis([-0.5 0.5])
set(gca,'FontSize',16)
% Plot average beta timecourse
beta_trace_avg = mean(beta_trace_sub,1);
beta_trace_std = std(beta_trace_sub,1);
lo = beta_trace_avg - beta_trace_std;
hi = beta_trace_avg + beta_trace_std;
% Define patch
subplot(1,2,2);
hp1 = patch([trl_time_induced, trl_time_induced(end:-1:1)], [lo, hi(end:-1:1)], 'r'); hold on;
h1 = plot(trl_time_induced,beta_trace_avg);
set(hp1, 'facecolor', [173 216 230]/255, 'edgecolor', 'none');
set(h1, 'color', [65 105 225]/255);
set(gcf,'color',[1 1 1])
set(gca,'FontSize',16)
xlabel('Time (s)'); ylabel('Beta power change from baseline')
xlim([min(trl_time_induced) max(trl_time_induced)])
mean_SNR = mean(beta_SNR_sub_MEGonly);
std_SNR = std(beta_SNR_sub_MEGonly);
mean_SNR_MEGonly = mean_SNR;
std_SNR_MEGonly = std_SNR;
beta_trace_avg_MEGonly = beta_trace_avg;
beta_trace_std_MEGonly = beta_trace_std;

% plot both EEG-MEG and MEGonly data together
figure;
lo_EEGMEG = beta_trace_avg_EEGMEG - beta_trace_std_EEGMEG;
hi_EEGMEG = beta_trace_avg_EEGMEG + beta_trace_std_EEGMEG;
hp2 = patch([trl_time_induced, trl_time_induced(end:-1:1)], [lo_EEGMEG, hi_EEGMEG(end:-1:1)], 'r'); hold on;
h2 = plot(trl_time_induced,beta_trace_avg_EEGMEG);
set(hp2, 'facecolor', [255 218 185]/255, 'edgecolor', 'none');
set(h2, 'color', [255 69 0]/255);
lo_MEGonly = beta_trace_avg_MEGonly - beta_trace_std_MEGonly;
hi_MEGonly = beta_trace_avg_MEGonly + beta_trace_std_MEGonly;
hp1 = patch([trl_time_induced, trl_time_induced(end:-1:1)], [lo_MEGonly, hi_MEGonly(end:-1:1)], 'r'); hold on;
h1 = plot(trl_time_induced,beta_trace_avg_MEGonly);
set(hp1, 'facecolor', [173 216 230]/255, 'edgecolor', 'none');
set(h1, 'color', [65 105 225]/255);
set(gcf,'color',[1 1 1])
set(gca,'FontSize',16)
xlabel('Time (s)'); ylabel('Beta power change from baseline')
xlim([min(trl_time_induced) max(trl_time_induced)])
alpha(0.4); % transparency value
legend({'MEG with EEG',' ','MEG only',' '},'FontSize',12)

% Plot SNR 
figure('Color','w'); 
bw = 0.8;
bar([mean_SNR_MEGonly,mean_SNR_EEGMEG],'BarWidth',bw); hold on
bar([0,mean_SNR_EEGMEG],'BarWidth',bw);
errorbar([mean_SNR_MEGonly,mean_SNR_EEGMEG],[std_SNR_MEGonly,std_SNR_EEGMEG],'k.')
set(gca,'FontSize',16); ylabel('Signal to Noise Ratio')
xticklabels({'MEG only','MEG with EEG'})
% Make scatters
scatter_vals_MEGonly = linspace(1-bw/2+0.1,1+bw/2-0.1,length(beta_SNR_sub_MEGonly));
s1 = scatter(scatter_vals_MEGonly,beta_SNR_sub_MEGonly);
s1.LineWidth = 0.6;
s1.MarkerEdgeColor = 'k';
s1.MarkerFaceColor = [1 93 182]/255;
scatter_vals_EEGMEG = linspace(2-bw/2+0.1,2+bw/2-0.1,length(beta_SNR_sub_EEGMEG));
s2 = scatter(scatter_vals_EEGMEG,beta_SNR_sub_EEGMEG);
s2.LineWidth = 0.6;
s2.MarkerEdgeColor = 'k';
s2.MarkerFaceColor = [200 95 20]/255;

% Get significance values for the difference between the group SNRs
[h,p] = ttest(beta_SNR_sub_MEGonly,beta_SNR_sub_EEGMEG)

% Number of trials
mean_Ntrials_EEGMEG = mean(Ntrials_sub_EEGMEG)
std_Ntrials_EEGMEG = std(Ntrials_sub_EEGMEG)
mean_Ntrials_MEGonly = mean(Ntrials_sub_MEGonly)
std_Ntrials_MEGonly = std(Ntrials_sub_MEGonly)

