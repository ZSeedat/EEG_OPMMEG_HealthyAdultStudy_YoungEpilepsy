% Function to convert a raw Cerca .cmeg file into fieldtrip format,
% for further analysis.
% Inputs are the raw MEG filename (including filepath).
% Output is fieldtrip data structure.
% Code edited from Nottingham scripts and Cerca Analysis app (Ryan Hill)
% for use at YE (Zelekha Seedat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Housekeeping
% % Add these paths in original script before calling function to save
% % computation time.
% 
% % Add fieldtrip path and load defaults
% addpath('C:\Documents\MATLAB\fieldtrip-20220212\fieldtrip-20220212')
% ft_defaults
% % Add path with data extraction code in it
% addpath('C:\Documents\MATLAB')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main function
function [FT_data_struct] = convert_2_fieldtrip(MEG_file)
Data = read_cMEG_data_split(MEG_file);
% Apply conversion factor to data
chan_info = tdfread([MEG_file(1:54),'channels.tsv'],'\t'); 
conv_fact = chan_info.nT0x2FV;
conv_fact(1:16,1:3) = repmat('1.0',16,1);
conv_fact = str2num(conv_fact);
Data.data = Data.data.*conv_fact;

% Get data sensor order
Data.label = cellstr(chan_info.name);

% Get trigger channels
Data.trigger = Data.data(1:16,:);

% Get position of sensors in helmet:
sens_pos = tdfread([MEG_file(1:54),'HelmConfig.tsv'],'\t'); 
sens_pos_id = sens_pos.Name;
idx = (1:64).*3-2;
sens_pos_order = sens_pos_id(idx,1:3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up data structs for X, Y and Z channels
% Get the indices for Y and Z channels
for i = 1:64 % 64 channels
    Z_inds(i) = 16+2*(i-1)+1;
    Y_inds(i) = 16+2*(i-1)+2;
end

% dataX.label = Data.label(X_inds);
dataY.label = Data.label(Y_inds);
dataZ.label = Data.label(Z_inds);
data_all.label = Data.label(16+1:end);

% [dataX.fsample,dataY.fsample,dataZ.fsample,data_all.fsample] = deal(Data.samp_frequency);
[dataY.fsample,dataZ.fsample,data_all.fsample] = deal(Data.samp_frequency);

% A single long trial which can be chopped into trials at a later stage
% dataX.trial{1} = Data.data(X_inds,:);
dataY.trial{1} = Data.data(Y_inds,:);
dataZ.trial{1} = Data.data(Z_inds,:);
data_all.trial{1} = Data.data(16+1:end,:);
% [dataX.time{1}, dataY.time{1}, dataZ.time{1}, data_all.time{1}] = deal([0:size(dataX.trial{1},2)-1]./dataX.fsample);
[dataY.time{1}, dataZ.time{1}, data_all.time{1}] = deal([0:size(dataY.trial{1},2)-1]./dataY.fsample);

% Move the channels into helmet order rather than rack order:
Helm_config = tdfread([MEG_file(1:54),'HelmConfig.tsv'],'\t');
helm_chan_order = cellstr(Helm_config.Sensor);
helm_chan_order = helm_chan_order(1:end-1); % the last entry isn't a sensor
% Remove extra channels because we only care about the order
chans = [1:length(helm_chan_order)/3].*3-2;
helm_chan_order = helm_chan_order(chans);
str_rack_label = char(dataY.label); % Y and Z channels are in the same order
inds_new = [];
for ch = 1:length(helm_chan_order)
    sens_to_find = helm_chan_order{ch};
    ch_ind = find(sum(str_rack_label(:,1:2) == sens_to_find(1:2),2)==2);
    inds_new = cat(1,inds_new,ch_ind);
end
dataY.label = dataY.label(inds_new);
dataZ.label = dataZ.label(inds_new);
data_all.label = [dataY.label;dataZ.label];
dataY.trial{1} = dataY.trial{1}(inds_new,:);
dataZ.trial{1} = dataZ.trial{1}(inds_new,:);
data_all.trial{1} = [dataY.trial{1};dataZ.trial{1}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare trigger events
% trigger_mat = Data.trigger;
% trig_time = Data.time;
% event = [];
% event_counter = 1;
% for trig_chan = 1:size(trigger_mat,1)
%     trig_str = ['trig',num2str(trig_chan)];
%     trig_n = trigger_mat(trig_chan,:);
%     diff_trig_n = diff(trig_n);
%     trig_up_inds = find(diff_trig_n > 3);
%     for i = 1:length(trig_up_inds)
%         event(event_counter).type = 'trigger';
%         event(event_counter).value = trig_str;
%         event(event_counter).sample = trig_up_inds(i)+1;
%         event(event_counter).duration = [];
%         event(event_counter).timestamp = trig_time(trig_up_inds(i));
%         event(event_counter).offset = 0;
%         event_counter = event_counter+1;
%     end
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare sensor locations
% pos_no = 1:64; % 64 sensor positions in helmet
% layout.pos = Helmet_info.lay.pos(pos_no,:);
% layout.width = Helmet_info.lay.width(pos_no,:);
% layout.height = Helmet_info.lay.height(pos_no,:);
% layout.label = Helmet_info.lay.label;
% 
% % For all channels - x, y, and z
% for chan_pos = 1:64 % 64 sensor positions in helmet
%     layout_all.pos(3*(chan_pos-1)+1,:) = Helmet_info.lay.pos(chan_pos,:);
%     layout_all.pos(3*(chan_pos-1)+2,:) = Helmet_info.lay.pos(chan_pos,:);
%     layout_all.pos(3*(chan_pos-1)+3,:) = Helmet_info.lay.pos(chan_pos,:);
% 
%     layout_all.width(3*(chan_pos-1)+1,:) = Helmet_info.lay.width(chan_pos,:);
%     layout_all.width(3*(chan_pos-1)+2,:) = Helmet_info.lay.width(chan_pos,:);
%     layout_all.width(3*(chan_pos-1)+3,:) = Helmet_info.lay.width(chan_pos,:);
% 
%     layout_all.height(3*(chan_pos-1)+1,:) = Helmet_info.lay.height(chan_pos,:);
%     layout_all.height(3*(chan_pos-1)+2,:) = Helmet_info.lay.height(chan_pos,:);
%     layout_all.height(3*(chan_pos-1)+3,:) = Helmet_info.lay.height(chan_pos,:);
% 
%     layout_all.label(3*(chan_pos-1)+1,1) = Helmet_info.lay.label(chan_pos);
%     layout_all.label(3*(chan_pos-1)+2,1) = Helmet_info.lay.label(chan_pos);
%     layout_all.label(3*(chan_pos-1)+3,1) = Helmet_info.lay.label(chan_pos);
% end
% 
% coreg_prmpt = questdlg('Do you have a coregistered sensor file to use?');
% if strcmpi('Yes',coreg_prmpt)
%     % Load transform here
%     maTrix = inputdlg({'R1','R2','R3',...
%         'R4','R5','R6',...
%         'R7','R8','R9',...
%         'T1','T2','T3'},...
%         'Transformation Matrix',[1 15],...
%         {'','','','','','','','','','','','',''});
% R = [str2num(maTrix{1}),str2num(maTrix{2}),str2num(maTrix{3});...
%     str2num(maTrix{4}),str2num(maTrix{5}),str2num(maTrix{6});...
%     str2num(maTrix{7}),str2num(maTrix{8}),str2num(maTrix{9})];
% T = [str2num(maTrix{10}),str2num(maTrix{11}),str2num(maTrix{12})];
% 
%     sens_info.pos = [R*Helmet_info.sens_pos']' + T;
%     sens_info.orsZ = [R*Helmet_info.sens_ors_Z']';
%     sens_info.orsY = [R*Helmet_info.sens_ors_Y']';
%     sens_info.orsX = [R*Helmet_info.sens_ors_X']';   
% else
%     sens_info.pos = Helmet_info.sens_pos;
%     sens_info.orsZ = Helmet_info.sens_ors_Z;
%     sens_info.orsY = Helmet_info.sens_ors_Y;
%     sens_info.orsX = Helmet_info.sens_ors_X;
% end
% % [gradX.coilpos,gradY.coilpos,gradZ.coilpos] = deal(sens_info.pos(pos_no,:).*100);
% [gradY.coilpos,gradZ.coilpos] = deal(sens_info.pos(pos_no,:).*100);
% % gradX.coilori = sens_info.orsX(pos_no,:);
% gradY.coilori = sens_info.orsY(pos_no,:);
% gradZ.coilori = sens_info.orsZ(pos_no,:);
% % gradX.label = dataX.label;
% gradY.label = dataY.label;
% gradZ.label = dataZ.label;
% % gradX.chanpos = gradX.coilpos;
% gradY.chanpos = gradY.coilpos;
% gradZ.chanpos = gradZ.coilpos;
% 
% % For all channels - both y and z
% for chan_pos = 1:64 % 64 sensor positions in helmet
% %     grad_all.coilpos(3*(chan_pos-1)+1,:) = Helmet_info.sens_pos(chan_pos,:); % No x-chans in our dual axis system
%     grad_all.coilpos(2*(chan_pos-1)+1,:) = Helmet_info.sens_pos(chan_pos,:);
%     grad_all.coilpos(2*(chan_pos-1)+2,:) = Helmet_info.sens_pos(chan_pos,:);
% 
% %     grad_all.coilori(3*(chan_pos-1)+1,:) = Helemt_info.sens_ors_X(chan_pos,:); % No x-chans in our dual axis system
%     grad_all.coilori(2*(chan_pos-1)+1,:) = Helmet_info.sens_ors_Y(chan_pos,:);
%     grad_all.coilori(2*(chan_pos-1)+2,:) = Helmet_info.sens_ors_Z(chan_pos,:);
% end
% grad_all.label = data_all.label;
% grad_all.chanpos = grad_all.coilpos;
% 
% % [gradX.units,gradY.units,gradZ.units,grad_all.units] = deal('cm');
% [gradY.units,gradZ.units,grad_all.units] = deal('cm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare Fieldtrip layout
% FT_data_struct_X.grad = gradX;
% FT_data_struct_X.data = dataX;
% FT_data_struct_X.layout = layout;

% FT_data_struct_Y.grad = gradY;
% FT_data_struct_Y.data = dataY;
% FT_data_struct_Y.layout = layout;

% FT_data_struct_Z.grad = gradZ;
% FT_data_struct_Z.data = dataZ;
% FT_data_struct_Z.layout = layout;

% FT_data_struct.grad = grad_all;
% FT_data_struct.data = data_all;
% FT_data_struct.layout = layout_all;

FT_data_struct.data_all = data_all;
FT_data_struct.dataZ = dataZ;
FT_data_struct.dataY = dataY;
FT_data_struct.trigger = Data.data(1:16,:);
FT_data_struct.sens_pos_order = sens_pos_order;
end

