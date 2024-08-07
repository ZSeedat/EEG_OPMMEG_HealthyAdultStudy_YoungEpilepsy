% Pre-process MEG in analysis app (Filter 1-150Hz, add the triggers in, notch,
% remove bad channels, mean field correct, remove bad trials, then export
% to workspace. Run convert_2_FT.m
% Sections for analysing both MEG and EEG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
addpath('R:\MATLAB\FieldTrip_analysis_scripts')
% Add fieldtrip path
addpath("R:\MATLAB\fieldtrip-20220212\fieldtrip-20220212")
ft_defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEG data
% Choose a MEG channel to look at:
for x = 1:length(data_fZ.label)
    match = strcmp('CR [Z]', data_fZ.label(x));
    if match == 1
        chan_ind = x;
    end
end

% Calculate PSD for 30s eyes open segment
MEG_Fs = data_fZ.fsample;
L = 5*MEG_Fs; % Length of segment
for t = 1:length(data_fZ.trial)
    fft_eyesOpen_megZ = fft(data_fZ.trial{1,1}(chan_ind,1*MEG_Fs:6*MEG_Fs));
    % This calculates the 2-sided spectrum
    P2_eo = abs(fft_eyesOpen_megZ/L);
    % but we only want the 1-sided:
    P1_eo = P2_eo(1:L/2+1);
    P1_eo(2:end-1) = 2*P1_eo(2:end-1);
    P1_eo_trials(t,:) = P1_eo;
    
    % Calculate PSD for 30s eyes closed segment
    fft_eyesClosed_megZ = fft(data_fZ.trial{1,1}(chan_ind,33*MEG_Fs+1:38*MEG_Fs));
    % This calculates the 2-sided spectrum
    P2_ec = abs(fft_eyesClosed_megZ/L);
    % but we only want the 1-sided:
    P1_ec = P2_ec(1:L/2+1);
    P1_ec(2:end-1) = 2*P1_ec(2:end-1);
    P1_ec_trials(t,:) = P1_ec;
end
% Average over trials:
mean_psd_eo = mean(P1_eo_trials,1);
% Smooth for plotting
eo_smooth = smoothdata(mean_psd_eo,'gaussian',10);
mean_psd_ec = mean(P1_ec_trials,1);
ec_smooth = smoothdata(mean_psd_ec,'gaussian',10);

% Plot
freqs = MEG_Fs*(0:(L/2))/L;
figure; plot(freqs,eo_smooth,'LineWidth',2); hold on;
plot(freqs,ec_smooth,'LineWidth',2);
xlim([0 30]); ylim([0 100])
xlabel('Frequency, Hz')
ylabel('P.S.D.')
title('MEG')
legend('Eyes open','Eyes closed')
set(gca,'FontSize',16)

