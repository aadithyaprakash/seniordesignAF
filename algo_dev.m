function [ total_weight ] = algo_dev( ecg_data, sampling_freq )
% We will be implementing the discrete dyadic wavelet transform using the
% methods from "A Wavelet-Based ECG Delineator: Evaluation on Standard
% Databases"
%
% If the prototype wavelet (mother wavelet) is the derivative of a
% smoothing function, then the Wavelet Transform of a signal at scale a is
% proportional to the derivative of the filtered signal w/ a smoothing
% impulse response @ scale a.
%
% Using literature, we have chosen to use a quadratic spline 

% %%% Dummy Data %%%%
% experiment_number = 0;
% experiment = ['ecg_test_', num2str(experiment_number), '.txt'];
% exp_data = load(experiment);
% t1 = 1:length(exp_data);

%%% Preprocessing ECG Data from MIT Database %%%%%%%%%%%%%%%%%%%%%%%%%
% 10800 points which corresponds to 30 s @ 360 sampling rate
[t1,exp_data] = rdsamp('mitdb/100', 1, 11000, 201); 

% % If we make this into a function that we call every 10s
% t1 = size(ecg_data);
% exp_data = ecg_data;

% Zero-mean
exp_data= exp_data - mean(exp_data(:));

% sampling frequency (what I assume the sampling frequency is)
f= sampling_freq;
% cutoff frequency
f_cutoff = 45; 
% cutoff frequency 2
f_cutoff_2 = .5;
% normalized cut off freq
fnorm =f_cutoff/(f/2); 
% normalized cut off freq 2
fnorm_2 =f_cutoff_2/(f/2); 

% First let's check the FFT to see what kind of frequencies we are dealing
% with
% fft_check(360, exp_data, 'FFT Original');

% Low pass + highpass filter combination to remove noise
% b = fir1(10, fnorm_2, 'high');
% exp_data_filtered = filter(b,1, exp_data);
% b = fir1(10, fnorm, 'low');
% exp_data_filtered = filter(b,1, exp_data);
% x = mean(exp_data_filtered)

%% REMOVE BASELINE WANDER USING MULTIRATE PROCESSING
% Decimate
% First low-pass filter at 45 Hz
b = fir1(10, fnorm, 'low');
lowpass_data = filter(b,1,exp_data);
% Then down-sample by 4
downsampled_data = downsample(lowpass_data,4);
% Then do a low pass filter at a low frequency
b = fir1(10, fnorm_2, 'low');
lowpass_data_2 = filter(b,1,downsampled_data);
% Interpolate
% upsample by 4
upsampled_data = upsample(lowpass_data_2,4);
% low pass filter 
wander = filter(b,1, upsampled_data);
% subtract this wander from the original signal
exp_data_filtered = exp_data - wander;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design FIR Filter as (1/8) * (d[n+2]+ 3d[n+1] + 3d[n-1]) where d is dirac
% impulse (NON-CAUSAL so first few ouputs will be garbage

a1 = 1;
b1 = [1/8, 3/8, 3/8, 1/8];

lowpass_filt_1 = filter(b1,a1,exp_data_filtered);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design FIR filter as 2 * (d[n+1] - d[n])

a2 = 1;
b2 = [0, 2, -2, 0];

highpass_filt = filter(b2,a2,exp_data);
highpass_filt_1 = filter(b2,a2,exp_data_filtered);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Interpolate signal by 2 before being fed into 2nd scale filters%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Low pass Butterworth filter of order 10 at Nyquist Frequency
[b_inter,a_inter] = butter(10,.98,'low'); 

% First thing to do is upsample by a factor of 2
y = upsample(lowpass_filt_1,2);

% Then run it through a low-pass filter at half the sampling freq
lowpass_filt_inter = filtfilt(b_inter,a_inter, y); 
% t2 = 1:length(lowpass_filt_inter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %2nd Scale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

highpass_filt_2 = filter(b2,a2,lowpass_filt_inter);
lowpass_filt_2 = filter(b1,a1, lowpass_filt_inter);
% 
% % To plot highpass_filt_2 output, I downsample by 2
% highpass_filt_2 = downsample(highpass_filt_2, 2);
% figure;
% plot(t1,exp_data_filtered,'-.',t1,highpass_filt_2,'-'), grid on
% legend('Original Data', 'Scale 2', 2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Interpolate signal by 2 before being fed into 3rd scale filters%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First thing to do is upsample by a factor of 2
y = upsample(lowpass_filt_2,2);

% Then run it through a low-pass filter at half the sampling freq
lowpass_filt_inter = filtfilt(b_inter,a_inter, y); 
% t2 = 1:length(lowpass_filt_inter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3rd Scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
highpass_filt_3 = filter(b2,a2,lowpass_filt_inter);
% lowpass_filt_3 = filter(b1,a1, lowpass_filt_inter);
highpass_filt_3 = downsample(highpass_filt_3, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pick out R-peaks
output = zeros(size(highpass_filt_1));
R_peaks = zeros(300,1); % Assume max heart rate is 50 beats every 10 seconds or 300 bpm.
R_index = 1;
R_counter = 0;
j = 20;

while (j < (length(highpass_filt_1)-10))
    if (highpass_filt_1(j-3) > .1 && highpass_filt_1(j+3) <-.1) && is_max(exp_data_filtered(j), exp_data_filtered(j+1), exp_data_filtered(j-1)) == 1 && exp_data_filtered(j) > 0.5
        output(j) = exp_data_filtered(j); % R-peak y position
        R_peaks(R_index) = j;
        j = j + 50; % Skip ahead to avoid double counting
        R_index = R_index + 1;    
        R_counter = R_counter + 1;
    else
        output(j) = 0;
        j = j + 1;
    end
end

%% Calculate variance of RR intervals
lengths = zeros(R_counter-1,1);
for j = 2:R_counter,
    lengths(j-1) = R_peaks(j)-R_peaks(j-1);
end

variance_RR = var(lengths);

%% Plot R-peaks
subplot(211);
plot(t1,exp_data_filtered,'-.',t1,highpass_filt_1,'-'), grid on
legend('Original Data', 'Scale 1', 2);
subplot(212);
plot(t1,exp_data_filtered,'-.',t1,output,'+','MarkerSize',5), grid on
legend('Data', 'R-Peaks', 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify P-waves
output_p = zeros(size(highpass_filt_1));
P_peaks = zeros(300,1); % Assume max heart rate is 50 beats every 10 seconds or 300 bpm.
P_index = 1;
P_counter = 0;
j = R_peaks(P_index) - 30;
counter = 0;

while(P_index < R_index)
    if j <= 30 || counter > 30
        output(j) = 0;
        counter = 0;
        if(P_index + 1 < R_index)
            P_index = P_index + 1;
            j = R_peaks(P_index) -30;
        else
            break;
        end
    elseif(highpass_filt_3(j-3) > 0 && highpass_filt_3(j+3) < 0) && (highpass_filt_1(j-3) < 0.3 && highpass_filt_1(j+3) > -0.3) && is_max(exp_data_filtered(j), exp_data_filtered(j+1), exp_data_filtered(j-1)) == 1
        if is_max(exp_data_filtered(j),exp_data_filtered(j+30), exp_data_filtered(j-30)) ==1
            output_p(j) = exp_data_filtered(j); % P-peak y position
            P_peaks(P_index) = j;
            if(P_index + 1 < R_index)
                P_index = P_index + 1;
                P_counter = P_counter + 1;
                j = R_peaks(P_index) -30;
                counter = 0;
            else
                P_counter = P_counter + 1;
                break;
            end
        else
            output(j) = 0;
            counter = counter + 1;
            j = j - 1; 
        end
    else
        output(j) = 0;
        counter = counter + 1;
        j = j - 1;
    end
end
figure;
subplot(211);
plot(t1,exp_data_filtered,'-.',t1,highpass_filt_3,'-'), grid on
legend('Original Data', 'Scale 3', 2);
subplot(212);
plot(t1,exp_data_filtered,'-.',t1,output_p,'+','MarkerSize',5), grid on
legend('Data', 'P-Peaks', 2);

%% Classification
% P_wave detection
if P_counter/R_counter > .90
    P_weight = 1;
elseif P_counter/R_counter > .60
    P_weight = .5;
else
    P_weight = 0;
end
% R-R interval variance
if variance_RR < 500
    RR_weight = 1;
elseif variance_RR < 1000
    RR_weight = .5;
else
    RR_weight = 0;
end
% Total weight from this section ??
total_weight = (.5 * P_weight) + (.5 * RR_weight);

