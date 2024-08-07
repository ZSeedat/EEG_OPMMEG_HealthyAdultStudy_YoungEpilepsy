% Script to read in EEG data for further analysis using FieldTrip
% Gets TFS and beta amplitude timecourse for highest SNR channel per
% subject and then averages over subjects
% Zelekha Seedat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
clear all; close all; clc;
% Add fieldtrip path
addpath("C:\Documents\MATLAB\fieldtrip-20220212\fieldtrip-20220212")
ft_defaults
% Choose directory to save results in
save_dir = uigetdir('M:\EEG MEG Project\Results','Choose directory to save results in');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in the data
% For all participants:
Pnum = {'P1','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12'}; % P2 has no EEG_MEG sync
for sub = 1:length(Pnum)
    clearvars -except Pnum sub save_dir
    [DataFile,DataPath]=uigetfile('.vhdr',['Select EEG data for ',Pnum{sub}],'M:\');
    % use Fieldtrip
    cfg = [];
    cfg.dataset = [DataPath,DataFile];
    cfg.demean = 'yes';
    cfg.detrend = 'yes';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49.5 50.5];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 150];
    cfg.reref = 'yes';
    cfg.refmethod = 'avg';
    cfg.refchannel = 'all';

    % Read in data
    data = ft_preprocessing(cfg);
    EEG_Fs = data.hdr.Fs;

    % At the minute all the data is one long trial
    data_all = data.trial{1};

    % Remove ECG channel
    ECG_chan = find(strcmp(data.label,'ECG'));
    data.label(ECG_chan) = [];
    data_all(ECG_chan,:) = [];
    Nchans = size(data_all,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import markers and create trials:
    % Extract MEG-EEG markers
    cfg.dataset = [DataPath,DataFile];
    event = ft_read_event(cfg.dataset,'detectflank',[]);

    % Write down the name of the EEG marker you would like to use. Look
    % for it in the 'event' variable in the MATLAB workspace.
    event_vals = {event.value};
    event_vals = char(event_vals{:});
    sel = ismember(event_vals,{'R128'}); % marker for stim cue
    trig_inds = find(sel);
    cue_marker_samp = [];
    for m = 1:length(trig_inds)
        cue_marker_samp(m) = event(trig_inds(m)).sample;
    end

    % User select the OptiTrack marker file:
    cd(DataPath)
    save_name = [DataFile(1:end-5),'_OT_trigs.mat'];
    load([DataPath,save_name]) % trig when finger abduction ends

    % Chop each trial
    offset = -1;
    seg_duration = 5; % 5s segments to average over
    data.trial = [];
    data.time = [];
    for t = 1:length(EEG_inds)
        data.trial{t} = data_all(:,EEG_inds(t)+EEG_Fs*offset:EEG_inds(t)+seg_duration*EEG_Fs+offset*EEG_Fs); % -10ms to 140ms around each marker
        data.time{t} = linspace(offset,seg_duration,size(data.trial{t},2));
        data.sampleinfo(t,:) = [EEG_inds(t)+EEG_Fs*offset, EEG_inds(t)+seg_duration*EEG_Fs+offset*EEG_Fs];
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
    EEG_inds = EEG_inds(good_trials);
    Ntrials = size(good_trials,2);
    disp([num2str(Ntrials) ' trials remain']) % Print to workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average over trials

    % % Calculate average for each trial
    % for ch = 1:size(data_all,1)
    %     for t = 1:length(data.trial)
    %         chan_trial = data.trial{1,t}(ch,:);
    %         chan_trial_bc(ch,t,:) = chan_trial-mean(chan_trial(end-1*EEG_Fs:end)); % baseline correct
    %     end
    %     mean_over_trials(ch,:) = mean(squeeze(chan_trial_bc(ch,:,:)),1);
    % end
    %
    % % Plot average for each trial
    % figure;
    % h=plot(data.time{1,1},mean_over_trials);
    % xlabel('Time, s'); ylabel('Trial average, uV'); title('All channels plotted')
    % legend(data.label,'NumColumns',3,'Location','northeastoutside')
    % for i = 1:numel(h)
    %     % Add a new row to DataTip showing the DisplayName of the line
    %     h(i).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Trace',repmat({h(i).DisplayName},size(h(i).XData)));
    % end
    % set(gca,'FontSize',16)
    %
    % % Plot all the trials for one channel to see if there is jitter in the
    % % stimulus artefact (jitter which could obscure the SEP)
    % chan = 'PO9';
    % chan_ind = find(strcmp(data.label,chan));
    % figure; plot(data.time{1,1},squeeze(chan_trial_bc(chan_ind,:,:)))
    % xlabel('Time, s'); ylabel('EEG, uV'); title(['Channel ',chan,': All trials plotted'])
    % set(gca,'FontSize',16)
    %
    % % Plot all the trials for C4 where we expect to see an SEP
    % chan = 'C4';
    % chan_ind = find(strcmp(data.label,chan));
    % figure; subplot(1,2,1); plot(data.time{1,1},squeeze(chan_trial_bc(chan_ind,:,:)))
    % xlabel('Time, s'); ylabel('EEG, uV'); title(['Channel ',chan,': All trials plotted'])
    % xlim([min(data.time{1,1}) max(data.time{1,1})]); set(gca,'FontSize',16)
    % subplot(1,2,2); plot(data.time{1,1},mean(squeeze(chan_trial_bc(chan_ind,:,:)),1))
    % xlabel('Time, s'); ylabel('EEG, uV'); title(['Channel ',chan,': Mean over trials'])
    % xlim([min(data.time{1,1}) max(data.time{1,1})]); set(gca,'FontSize',16)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Use Filedtrip to average as a sanity check

    % cfg = [];
    % tlk = ft_timelockanalysis(cfg,data);
    %
    % % Use EEG layout to plot results
    % cd('C:\Users\MEG Group\OneDrive - Young Epilepsy\Documents\MATLAB\fieldtrip-20220212\fieldtrip-20220212\template\layout')
    % layout = load('layout_Nottscap64.mat'); % 64 chan arrangement, used in BrainCap64

    % cfg = [];
    % cfg.layout = layout.lay;
    % figure; ft_plot_layout(layout.lay)
    % figure; ft_multiplotER(cfg,tlk);
    % cfg.xlim = [0.02 0.03];
    % cfg.interactive = 'yes';
    % figure; ft_topoplotER(cfg,tlk); colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at time-frequency analysis for induced response
    % % Use sliding window analysis
    % cfg              = [];
    % cfg.output       = 'pow';
    % cfg.channel      = 'all';
    % cfg.method       = 'mtmconvol';
    % cfg.taper        = 'hanning';
    % cfg.foi          = 4:2:40;                         % analysis 2 to 40 Hz in steps of 2 Hz
    % cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    % cfg.toi          = 0:0.05:5;                    % the time window "slides" from -0.5 to 1.5 in 0.05 sec steps
    % cfg.pad = 'nextpow2';
    % TFRhann_MNS = ft_freqanalysis(cfg,data);    % visual stimuli
    %
    % % Plot
    % cfg = [];
    % cfg.baseline     = [4 5];
    % cfg.baselinetype = 'absolute';
    % cfg.showlabels   = 'yes';
    % cfg.layout       = layout.lay;
    % cfg.xlim = 'maxmin'; cfg.ylim = 'maxmin'; cfg.zlim = 'maxmin';
    % figure; ft_multiplotTFR(cfg, TFRhann_MNS);
    %
    % % Wavelet analysis (morlet)
    % cfg = [];
    % cfg.channel    = 'all';
    % cfg.method     = 'wavelet';
    % cfg.width      = 7;
    % cfg.output     = 'pow';
    % cfg.foi        = 1:2:40;
    % cfg.toi        = 0:0.05:5;
    % cfg.fsample = data.fsample;
    % TFRwave_MNS = ft_freqanalysis(cfg, data);
    %
    % % Plot
    % cfg = [];
    % cfg.baseline     = [4 5];
    % cfg.baselinetype = 'absolute';
    % cfg.marker       = 'on';
    % cfg.layout       = layout.lay;
    % cfg.xlim = 'maxmin'; cfg.ylim = 'maxmin'; cfg.zlim = 'maxmin';
    % figure; ft_multiplotTFR(cfg, TFRwave_MNS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot time-frequency spectrogram
    highpass = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 65 70];% 75 80 85 90 95 100 105 110];
    lowpass = [4 6 8 10 13 20 25 30 35 40 45 50 55 60 65 70 75 80];% 85 90 95 100 105 110 115 120];
    fre = highpass + ((lowpass - highpass)./2);
    % Control window
    % duration = input('Please type in the trial duration in seconds, then press enter.');
    duration = 5; % 5 second trials
    MRBD_win = ([-1, 0]-offset).*EEG_Fs+1; % desync window is 1 to 0s before movement cessation
    con_win = [duration-0.5 duration].*EEG_Fs;
    trl_time_induced = linspace(offset,duration+offset,duration*EEG_Fs);

    % Initialise figure for plotting
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
    % Use EEG layout to plot results
    cd('C:\Users\MEG Group\OneDrive - Young Epilepsy\Documents\MATLAB\fieldtrip-20220212\fieldtrip-20220212\template\layout')
    layout = load('layout_Nottscap64.mat'); % 64 chan arrangement, used in BrainCap64

    % Loop over each sensor
    TFS = [];
    data_all_labels = data.cfg.previous.channel;
    data_all_labels(ECG_chan) = [];
    Nchans = length(data_all_labels);
    for n = 1:Nchans
        chan_idx = find(startsWith(data_all_labels,layout.lay.label{n}));
        if chan_idx
            EEG_ch = data_all(chan_idx,:)';
            % Filter data within bands and calculate envelope
            EEG_ch_fb = zeros(length(EEG_ch),length(fre));
            for fb = 1:length(highpass)
                fprintf('\n Band %u/%u ',fb,length(highpass))
                filt_dat = nut_filter3_nottsapp(EEG_ch,'butter','bp',3,highpass(fb),lowpass(fb),EEG_Fs,1);
                EEG_ch_fb(:,fb) = abs(hilbert(filt_dat));
            end
            EEG_mean = zeros(duration*EEG_Fs,length(fre));
            % Chop data into trials
            for fb = 1:length(highpass)
                EEG_ch_fb_trials = [];
                for i = 1:Ntrials
                    EEG_ch_fb_trials = cat(1,EEG_ch_fb_trials,EEG_ch_fb(EEG_inds(i)+offset*EEG_Fs+1:EEG_inds(i)+offset*EEG_Fs+(duration*EEG_Fs),fb));
                end
                EEG_ch_filt = reshape(EEG_ch_fb_trials,duration*EEG_Fs,Ntrials);
                % Average across trials
                EEG_mean(:,fb) = mean(EEG_ch_filt,2);
            end
            meanrest = mean(EEG_mean(con_win(1):con_win(2),:),1);
            meanrestmat = repmat(meanrest,size(EEG_mean,1),1);
            TFS = (EEG_mean'-meanrestmat')./meanrestmat';

            % plot
            figure(200)
            if exist('layout')
                p = uipanel('Title',layout.lay.label(n),...
                    'Units','normalized','Position',...
                    [-0.125+layout.lay.pos(n,1)/580 (-0.075+layout.lay.pos(n,2)/480)*0.95 0.07 0.09],...
                    'BackgroundColor','w','BorderType','none');
                ax = axes(p);

            else
                subplot(ceil(sqrt(Nchans)),ceil(sqrt(Nchans)),n);
            end
            pcolor(trl_time_induced,fre,TFS);shading interp
            axis fill
            caxis([-0.5 0.5])
            set(gcf,'color',[1 1 1])
            if ~exist('layout')
                xlabel('Time (s)');ylabel('Frequency (Hz)')
                title(layout.lay.label(n))
            end
            axis off
            title(sprintf('Ch %s',num2str(n)))
            pause(0.1)
        else
        end
    end
    p = uipanel('Title','Axes',...
        'Units','normalized','Position',...
        [-0.125+layout.lay.pos(n,1)/580-0.5 (-0.075+layout.lay.pos(n,2)/480)*0.95 0.07 0.09],...
        'BackgroundColor','w','BorderType','none');
    ax = axes(p); ylim([min(fre) max(fre)])
    xlim([min(trl_time_induced) max(trl_time_induced)])
    xlabel('Time (s)'); ylabel('Freq (Hz)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select the sensor with the greatest beta-band SNR to save out
    beta_highpass = 13;
    beta_lowpass = 30;
    PMBR_win = ([1, 2]-offset).*EEG_Fs; % rebound window is 1 to 2s after movement cessation
    MRBD_win = ([-1, 0]-offset).*EEG_Fs+1; % desync window is 1 to 0s before movement cessation
    con_win = [duration-0.5 duration].*EEG_Fs;
    for n = 1:Nchans
        chan_idx = find(startsWith(data_all_labels,layout.lay.label{n}));
        if chan_idx
            EEG_ch = data_all(chan_idx,:)';

            % Filter data into beta bands and calculate envelope
            filt_dat = nut_filter3_nottsapp(EEG_ch,'butter','bp',3,beta_highpass,beta_lowpass,EEG_Fs,1);
            EEG_ch_beta = abs(hilbert(filt_dat));

            % Chop data into trials
            EEG_ch_beta_trials = [];
            for i = 1:Ntrials
                EEG_ch_beta_trials = cat(1,EEG_ch_beta_trials,EEG_ch_beta(EEG_inds(i)+offset*EEG_Fs+1:EEG_inds(i)+offset*EEG_Fs+(duration*EEG_Fs)));
            end
            EEG_ch_beta_filt = reshape(EEG_ch_beta_trials,duration*EEG_Fs,Ntrials);
            % Average across trials
            EEG_mean = mean(EEG_ch_beta_filt,2);

            meanrest = mean(EEG_mean(con_win(1):con_win(2)));
            meanrestmat = repmat(meanrest,size(EEG_mean,1),1);
            beta_trace(n,:) = (EEG_mean'-meanrestmat')./meanrestmat';
            beta_chan_SNR(n) = (mean(beta_trace(n,PMBR_win(1):PMBR_win(2)))-mean(beta_trace(n,MRBD_win(1):MRBD_win(2))))/std(beta_trace(n,MRBD_win(1):MRBD_win(2)));

            % plot
            figure(300)
            if exist('layout')
                p = uipanel('Title',layout.lay.label(n),...
                    'Units','normalized','Position',...
                    [-0.125+layout.lay.pos(n,1)/580 (-0.075+layout.lay.pos(n,2)/480)*0.95 0.07 0.09],...
                    'BackgroundColor','w','BorderType','none');
                ax = axes(p);

            else
                subplot(ceil(sqrt(Nchans)),ceil(sqrt(Nchans)),n);
            end
            plot(trl_time_induced,beta_trace(n,:));
            axis fill
            set(gcf,'color',[1 1 1])
            if ~exist('layout')
                xlabel('Time (s)');ylabel('Frequency (Hz)')
                title(layout.lay.label(n))
            end
            axis off
            title(sprintf('Ch %s',num2str(n)))
            ylim([-1 1])
            pause(0.1)
        else
        end
    end
    % Sensor with max SNR in the left centro-parietal electrodes
    relevant_electrodes = [26, 27, 28, 36, 37, 38, 46, 47, 48];
    [max_SNR,max_SNR_idx] = max(beta_chan_SNR(relevant_electrodes));
    chan = data_all_labels{relevant_electrodes(max_SNR_idx)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot TFS and beta timecourse for selected sensor
    chan_idx = find(strcmp(chan,data_all_labels));
    MRBD_win = ([-1, 0]-offset).*EEG_Fs+1; % desync window is 1 to 0s before movement cessation
    con_win = [duration-0.5 duration].*EEG_Fs;
    offset = -1;
    EEG_ch = data_all(chan_idx,:)';
    % Filter data within bands and calculate envelope
    EEG_ch_fb = zeros(length(EEG_ch),length(fre));
    for fb = 1:length(highpass)
        fprintf('\n Band %u/%u ',fb,length(highpass))
        filt_dat = nut_filter3_nottsapp(EEG_ch,'butter','bp',3,highpass(fb),lowpass(fb),EEG_Fs,1);
        EEG_ch_fb(:,fb) = abs(hilbert(filt_dat));
    end
    EEG_mean = zeros(duration*EEG_Fs,length(fre));
    % Chop data into trials
    for fb = 1:length(highpass)
        EEG_ch_fb_trials = [];
        for i = 1:Ntrials
            EEG_ch_fb_trials = cat(1,EEG_ch_fb_trials,EEG_ch_fb(EEG_inds(i)+offset*EEG_Fs+1:EEG_inds(i)+offset*EEG_Fs+(duration*EEG_Fs),fb));
        end
        EEG_ch_filt = reshape(EEG_ch_fb_trials,duration*EEG_Fs,Ntrials);
        % Average across trials
        EEG_mean(:,fb) = mean(EEG_ch_filt,2);
    end
    meanrest = mean(EEG_mean(con_win(1):con_win(2),:),1);
    meanrestmat = repmat(meanrest,size(EEG_mean,1),1);
    TFS = (EEG_mean'-meanrestmat')./meanrestmat';
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
    filt_dat = nut_filter3_nottsapp(EEG_ch,'butter','bp',3,beta_highpass,beta_lowpass,EEG_Fs,1);
    EEG_ch_beta = abs(hilbert(filt_dat));

    % Chop data into trials
    EEG_ch_beta_trials = [];
    for i = 1:Ntrials
        EEG_ch_beta_trials = cat(1,EEG_ch_beta_trials,EEG_ch_beta(EEG_inds(i)+offset*EEG_Fs+1:EEG_inds(i)+offset*EEG_Fs+(duration*EEG_Fs)));
    end
    EEG_ch_beta_filt = reshape(EEG_ch_beta_trials,duration*EEG_Fs,Ntrials);
    % Average across trials
    EEG_mean = mean(EEG_ch_beta_filt,2);

    meanrest = mean(EEG_mean(con_win(1):con_win(2)));
    meanrestmat = repmat(meanrest,size(EEG_mean,1),1);
    beta_trace = (EEG_mean'-meanrestmat')./meanrestmat';
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
    cd(save_dir)
    mkdir([save_dir,'\',Pnum{sub}])
    cd([save_dir,'\',Pnum{sub}])
    save('TFS.mat','TFS','trl_time_induced','fre','chan','Ntrials')
    save('beta_trace.mat','beta_trace','trl_time_induced','chan','Ntrials')
    save('beta_chan_SNR.mat','beta_chan_SNR','chan','Ntrials')

    drawnow
    pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot group average result:
% For all participants:
Pnum = {'P1','P2','P3','P5','P6','P7','P10','P11','P12'}; % P4 has no EEG_MEG sync, P8 and P9 had no rebound
save_dir = 'R:\EEG MEG Project\Results\EEGMEG_EEGdata\R_index_abduction';
clearvars -except save_dir Pnum
for sub = 1:length(Pnum)
    cd([save_dir,'\',Pnum{sub}])
    load("TFS.mat")
    TFS_sub(:,:,sub) = TFS;
    load("beta_trace.mat")
    beta_trace_sub(sub,:) = beta_trace;
    load("beta_chan_SNR.mat")
    beta_SNR_sub(sub) = beta_chan_SNR;
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
mean_SNR = mean(beta_SNR_sub);
std_SNR = std(beta_SNR_sub);
mean_SNR_EEGMEG = mean_SNR;
std_SNR_EEGMEG = std_SNR;
beta_trace_avg_EEGMEG = beta_trace_avg;
beta_trace_std_EEGMEG = beta_trace_std;

