%% Main script for the estimation of running temporal parameters
% Author: Rachele Rossanigo, r.rossanigo@phd.uniss.it
% May 2024

%% ----------------------------Data load-----------------------------------
% Load templates
load temple_library.mat
% template_library contains three speed-specific templates (at low,
% moderate, and fast running speeds). Each template contains two Nx3
% matrixes (one for the accelerations in m/s^2 and one for the angular
% velocities in rad/s): first column corresponds to the anteroposterior 
% axis, second column to the mediolateral axis (towards subject's right),
% third column to the vertical axis. Each template contains its initial 
% contact (IC) in frames, final contact (FC) in frames, and its sampling
% frequency. 

% Upload data to analyze
...
Acc = ...; % accelearion (m/s^2)
Gyr = ...; % angular velocity (deg/s)
fs=...; % sampling frequency (Hz)

% Acc and Gyr must include anteroposterior axis data in first column,
% mediolateral axis (towards subject's right) data in second column, and
% vertical axis data in third column. 

% Sampling frequency
fs = ...

% Preprocessing
% Removal of the gyroscope bias estimated in a static acquisition
...

% Identification of the running interval to analyze
runningPortion = ...; % interval in samples

%% ------------------Identification of mid-swing instants------------------
runningPortion = round(runningPortion);
Gyr_ML = Gyr(runningPortion,2); % Selection of the mediolateral angular velocity (for instance in this case the second column)

% Power spectral density of the mediolateral angular velocity
window = round(fs/0.5); noverlap = window/2; NFFT = fs/0.25; % Set parameters for the power spectral density 
[pxx,f] = pwelch(Gyr_ML-mean(Gyr_ML),hamming(window),noverlap,NFFT,fs); 
figure,plot(f,pxx), title('Power Spectral Density of mediolateral angular velocity')
ylabel('(rad/s^2)/Hz'), xlabel('Hz')

% Definition of the parametric cutoff frequency based on the stride frequency
[~, in_f_stride] = max(pxx);
f_stride = f(in_f_stride);
cutoff = 1.5*f_stride;

% Filtering to enhance the angular velocity peaks corresponding to mid-swing instants
Wn = cutoff/(fs/2);
[B,A] = butter(2, Wn, 'low'); % Butterworth 2nd-order low-pass filter
Gyr_filt=filtfilt(B,A, Gyr_ML); % filtered signal with enhanced peaks

% Identification of the mid-swings instants
[pp, ms] = findpeaks(Gyr_filt, 'MinPeakHeight', 0, 'MinPeakDistance', 0.3*fs);
figure, plot(Gyr_filt), hold on, plot(ms, pp,'*'), 
xlabel('Frames (#)'), ylabel('Angular velocity (rad/s)')
legend('Filtered mediolateral angular velocity', 'Mid-swing instants')

%% ----------------------------Pre-processing------------------------------ 
% low-pass filtering to remove high-frequency noise
filterOrder = 4; cutoff = 35; Wn = cutoff/(fs/2); type = 'low';
[b,a] = butter(filterOrder,Wn,type);  % Butterworth 4th-order low-pass filter
Acc = filtfilt(b,a,Acc);
Gyr = filtfilt(b,a,Gyr);

%% -------------------------IC and FC detection----------------------------
IC = zeros(length(ms)-1,2); FC =  zeros(length(ms)-1,2); % initialization 
% IC and FC are Mx2 matrix, where M is the number of mid-swing instants
% minus one. The first column will contain the frames of identified
% initial and final contact, while the second column will contain the ID
% of the selected template (i.e., 1, 2, or 3). 

% Mid-swing to mid-swing segmentation
for k=1:length(ms)-1
    cycle_window = ms(k):ms(k+1);
    if length(cycle_window)>= 0.3*fs && length(cycle_window)<= 1*fs % check on mid-swing to mid-swing distance: a running gait cycle should last between 0.3 s and 1 s 
        % Template-based estimation of initial and final contact instants
        [IC(k,:), FC(k,:)]=template_based_IC_FC_detection(runningPortion, cycle_window, Acc, Gyr, template_library, fs);
    else
        FC(k,:) = [NaN, NaN];
        IC(k,:) = [NaN, NaN];
    end
end

% Plot
figure, p1 = plot(Gyr(:,2)); p2=xline(IC(:,1),'--g'); p3=xline(FC(:,1),'--m');                         
legend([p1,p2(1),p3(1)],{'Mediolateral angular velocity', 'Initial contact','Final contact'}),
ylabel('Angular velocity (rad/s)'), xlabel('Frame (#)')

%% -------------------Computation of temporal parameters-------------------
IC_t = IC(:,1)/fs; % IC in s
FC_t = FC(:,1)/fs; % FC in s
stride = diff(IC_t(:,1)); % running gait cycle duration or stride duration
stance = FC_t(:,1)-IC_t(:,1); % stance duration (foot in contact with the ground)
swing = IC_t(2:end,1)-FC_t(1:end-1,1); % swing duration (foot not in contact with the ground)
duty_factor = stance(1:end-1)./stride*100; % duty factor


                       


