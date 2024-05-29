%% TEMPLATE-BASED METHOD FOR THE DETECTION OF INITIAL AND FINAL CONTACT THROUGH DYNAMIC TIME WARPING (DTW)

function [IC, FC] = template_based_IC_FC_detection(runningPortion, cycle_window, Acc, Gyr, templates, fs_data)

% INPUT:
% runningPortion: complete running interval to analyze (samples)
% cycle_window: running gait cycle interval between two consecutive mid-swing instants (samples)
% Acc: acceleration signals m/s^2 (x,y,z), with x -> anteroposterior axis, y -> mediolateral axis (towards subject'sright), z -> vertical axis
% Gyr: gyroscope signals rad/s (x,y,z), with x -> anteroposterior axis, y -> mediolateral axis (towards subject'sright), z -> vertical axis
% fs_data: sampling frequency of the dataset to analyze
% templates: structure containing three annotated templates

% OUTPUT:
% IC and FC:
% IC(1) and FC(1) contain the initial and final contact in sample of the current running gait to analyze
% IC(2) and FC(2) contain the id of the machted template (1, 2, or 3)

% Author: Rachele Rossanigo, r.rossanigo@phd.uniss.it

%% ---------------------------IC DETECTION--------------------------------- 
% Computation of acceleration norm
Acc_norm = sqrt(Acc(:,1).^2+Acc(:,2).^2+Acc(:,3).^2);

% Segmentatation within the current running gait cycle
acc_norm_signal=Acc_norm(runningPortion(cycle_window),:);

templates_ids = string(fieldnames(templates)); % templates names
N = length(templates_ids); % templates number
dtw_distances = zeros(N,1); IC = []; % initialization

%% DTW-based comparison with the templates
for i = 1:N
    % i-th acceleration template
    template=char(templates_ids{i});
    acc_template = templates.(template).Acc;
    fs_template = templates.(template).fs;

    % Filtering
    filterOrder = 4;
    cutoff = 35;
    Wn = cutoff/(fs_template/2);
    type = 'low';
    [b,a] = butter(filterOrder,Wn,type);
    acc_template = filtfilt(b,a,acc_template);

    % Template acceleration norm
    acc_norm_template = sqrt(acc_template(:,1).^2+acc_template(:,2).^2+acc_template(:,3).^2);

    % DTW-based similarity between the running gait cycle and the i-th template
    [dist,ix,iy] = dtw(acc_norm_template, acc_norm_signal, 'euclidean');
    % dtw Matlab function stretches two vectors, acc_norm_template and 
    % acc_norm_signal, onto a common set of instants such that dist (i.e., 
    % the sum of the Euclidean distances between corresponding points)
    % is smallest. ix and iy returns the warping path (set of instants)
    % such that acc_norm_template(ix) and acc_norm_signal(iy) have the 
    % smallest possible dist between them.
    
    dtw_distances(i)=dist; % dtw-based similarity measure (Euclidean distance) between the i-th template and the running gait cycle 
    dtw_distortions.(template).ix=ix; % dtw-based distortion of the template
    dtw_distortions.(template).iy=iy; % dtw-based distotion of the running gait cycle to analyze
end

%% Template selection
[~,idx]=min(dtw_distances); % the most similar template corresponds to the minimum distance
closestMatch=char(templates_ids{idx}); % selection of the most similar template through the id

%% Computation of IC of the analyzed running cycle
% IC of the selected annotated template is known
IC_template =templates.(closestMatch).IC; 

% Define the sample of the warped signals corresponding to IC
IC_dtw_vett=find(dtw_distortions.(closestMatch).ix==int64(IC_template)); % find the sample(s) of the warped signals corresponding to IC 
if length(IC_dtw_vett)==1
    IC_dtw = IC_dtw_vett;
else
    IC_dtw=IC_dtw_vett(round(length(IC_dtw_vett)/2));
end

% Find the sample of the original running gait cycle corresponding to IC
IC_signal = dtw_distortions.(closestMatch).iy(IC_dtw); 
IC(1)=runningPortion(cycle_window(1))+ IC_signal; % translate the sample into the trial time vector
IC(2)=idx; % save the id of the selected template

%% Plot
acc_template_match = templates.(closestMatch).Acc; 
ix = dtw_distortions.(closestMatch).ix; 
iy = dtw_distortions.(closestMatch).iy; 
fs_template = templates.(closestMatch).fs;

acc_norm_template_match = sqrt(acc_template_match(:,1).^2+acc_template_match(:,2).^2+acc_template_match(:,3).^2);
    
filterOrder = 4;
cutoff = 35;
Wn = cutoff/(fs_template/2);
type = 'low';
[b,a] = butter(filterOrder,Wn,type);
acc_norm_template_match = filtfilt(b,a,acc_norm_template_match);

time_template = 0:1/fs_template:(length(acc_norm_template_match)-1)/fs_template; 
time_signal = 0:1/fs_data:(length(acc_norm_signal)-1)/fs_data; 

figure, plot(acc_norm_template_match), hold on, plot(acc_norm_signal), xline(IC_template),
title ('Acceleration norms in their own original time intervals'), legend('template','signal')
ylabel('m/s^2'), xlabel('Frame (#)')

figure, plot(acc_norm_template_match(ix)), hold on,  plot(acc_norm_signal(iy)), xline(IC_dtw), 
title ('DTW-aligned acceleration norms'), legend('template','signal')
ylabel('m/s^2'), xlabel('Frame (#)')

figure, plot(time_template, acc_norm_template_match), hold on, plot(time_signal, acc_norm_signal),
xline(time_template(IC_template),'b'), xline(time_signal(IC_signal),'r'),
title ('Acceleration norms in their own original time intervals with their ICs'), legend('template','signal')
ylabel('m/s^2'), xlabel('Time (s)')

%% ----------------------------FC DETECTION--------------------------------
% Segmentation of mediolateral angular velocity of the current running gait cycle
gyr_ml_signal=Gyr(runningPortion(cycle_window),2); 

dtw_distances=zeros(N,1); FC = []; % initialization

%% DTW-based comparison with the templates
for i = 1:N
    % i-th angular velocity template
    template=char(templates_ids{i});
    gyr_template = templates.(template).Gyr;
    fs_template = templates.(template).fs;

    % Filtering
    filterOrder = 4;
    cutoff = 35;
    Wn = cutoff/(fs_template/2);
    type = 'low';
    [b,a] = butter(filterOrder,Wn,type);
    gyr_template = filtfilt(b,a,gyr_template);

    gyr_ml_template = gyr_template(:,2); % mediolateral angular rate of the template

    % DTW-based similarity between the running gait cycle and the i-th template
    [dist,ix,iy] = dtw(gyr_ml_template, gyr_ml_signal, 'euclidean'); 
    % dtw Matlab function stretches two vectors, gyr_ml_template and 
    % gyr_ml_signal, onto a common set of instants such that dist (i.e., 
    % the sum of the Euclidean distances between corresponding points)
    % is smallest. ix and iy returns the warping path (set of instants)
    % such that gyr_ml_template(ix) and gyr_ml_signal(iy) have the smallest 
    % possible dist between them.

    dtw_distances(i)=dist; % dtw-based similarity (Euclidean distance) between the i-th template and the current running gait cycle 
    dtw_distortions.(template).ix=ix; % dtw-based distortion of the template
    dtw_distortions.(template).iy=iy; % dtw-based distotion of the running cycle to analyze
end

%% Template selection
[~,idx]=min(dtw_distances); % the most similar template corresponds to the minimum distance
closestMatch=char(templates_ids{idx}); % selection of the most similar template through the id

%% Computation of FC of the analyzed running cycle
% FC of the selected annotated template is known
FC_template =templates.(closestMatch).FC; 

% Define the sample of the warped signals corresponding to FC
FC_dtw_vett=find(dtw_distortions.(closestMatch).ix==int64(FC_template)); % find the sample(s) of the warped signals corresponding to FC 
if length(FC_dtw_vett)==1
    FC_dtw = FC_dtw_vett;
else
    FC_dtw=FC_dtw_vett(round(length(FC_dtw_vett)/2)); % this definition of FC_dtw can be modified in the optional check below
end

% Find the sample of the original running gait cycle corresponding to IC
FC_signal = dtw_distortions.(closestMatch).iy(FC_dtw); 

%--------------------------------------------------------------------------
% Optional: if there are more samples in the warped time vector
% potentially correspoding to FC, a further check can be implemented 
% assuming a minimum value of accepatable stance duration (i.e., 0.10 s)
iii=round(length(FC_dtw_vett)/2);
while (FC_signal-IC_signal)/fs_data<0.1 && iii < length(FC_dtw_vett)+1
    FC_dtw=FC_dtw_vett(iii);
    FC_signal = dtw_distortions.(closestMatch).iy(FC_dtw);
    iii=iii+1;
end
%--------------------------------------------------------------------------

FC(1)=runningPortion(cycle_window(1))+ FC_signal; % translate the sample into the trial time vector
FC(2)=idx; % save the template matching info

%% Plot
gyr_ML_template_match = templates.(closestMatch).Gyr(:,2);
ix = dtw_distortions.(closestMatch).ix; 
iy = dtw_distortions.(closestMatch).iy; 
fs_template = templates.(closestMatch).fs;

filterOrder = 4;
cutoff = 35;
Wn = cutoff/(fs_template);
type = 'low';
[b,a] = butter(filterOrder,Wn,type);
gyr_ML_template_match = filtfilt(b,a,gyr_ML_template_match);

time_template = 0:1/fs_template:(length(gyr_ML_template_match)-1)/fs_template; 
time_signal = 0:1/fs_data:(length(gyr_ml_signal)-1)/fs_data; 

figure, plot(gyr_ML_template_match), hold on, plot(gyr_ml_signal), xline(FC_template),
title ('Mediolateral angular velocities in their own original time intervals'), legend('template','signal')
ylabel('rad/s'), xlabel('Frame (#)')

figure, plot(gyr_ML_template_match(ix)), hold on,  plot(gyr_ml_signal(iy)), xline(FC_dtw), 
title ('DTW-aligned mediolateral angular velocities'), legend('template','signal')
ylabel('rad/s'), xlabel('Frame (#)')

figure, plot(time_template, gyr_ML_template_match), hold on, plot(time_signal,gyr_ml_signal),
xline(time_template(FC_template),'b'), xline(time_signal(FC_signal),'r'),
title ('Mediolateral angular velocities in their own original time intervals with their ICs'), legend('template','signal')
ylabel('rad/s'), xlabel('Time (s)')


end
