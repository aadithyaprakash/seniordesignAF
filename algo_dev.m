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
% 3600 points which corresponds to 10 s @ 360 sampling rate
[t1,exp_data] = rdsamp('mitdb/100', 1,1800,200); 

% Zero-mean
exp_data = exp_data - mean(exp_data(:));

% sampling frequency (what I assume the sampling frequency is)
f=360;
% cutoff frequency
f_cutoff = 40; 
% normalized cut off freq
fnorm =f_cutoff/(f/2); 

% Low pass Butterworth filter of order 10
[b_but,a_but] = butter(10,fnorm,'low'); 
% filtering the raw data w/ low-pass filter w/ cutoff at 40 Hz
exp_data_filtered = filtfilt(b_but,a_but,exp_data); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design FIR Filter as (1/8) * (d[n+2]+ 3d[n+1] + 3d[n-1]) where d is dirac
% impulse (NON-CAUSAL so first few ouputs will be garbage

a1 = 1;
b1 = [1/8, 3/8, 3/8, 1/8];

lowpass_filt_1 = filter(b1,a1,exp_data_filtered);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design FIR filter as 2 * (d[n+1] - d[n])

a2 = 1;
b2 = [0, 2, -2, 0];

highpass_filt = filter(b2,a2,exp_data);
highpass_filt_1 = filter(b2,a2,exp_data_filtered);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick out R-peaks
output = zeros(size(highpass_filt_1));
R_peaks = zeros(50:1); % Assume max heart rate is 50 beats every 10 seconds or 300 bpm.
R_index = 1;
j = 2;

while (j < (length(highpass_filt_1)-1))
    if (highpass_filt_1(j-1) > 0 && highpass_filt_1(j+1) <0) && (exp_data_filtered(j)) > .5,
        output(j) = 1; % R-peak
        R_peaks(R_index) = j;
        j = j + 50; % Skip ahead to avoid double counting
        R_index = R_index + 1;
    else
        output(j) = 0;
        j = j + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First Scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(211);
plot(t1,exp_data_filtered,'-.',t1,highpass_filt_1,'-'), grid on
legend('Data', 'Scale-1 Data', 2);
subplot(212);
plot(t1,exp_data_filtered,'-.',t1,output,'-'), grid on
legend('Data', 'Peaks', 2);

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
highpass_filt_2 = downsample(highpass_filt_2, 2);
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
lowpass_filt_3 = filter(b1,a1, lowpass_filt_inter);
highpass_filt_3 = downsample(highpass_filt_3, 4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify P-waves
output_p = zeros(size(highpass_filt_1));
P_peaks = zeros(50:1); % Assume max heart rate is 50 beats every 10 seconds or 300 bpm.
P_index = 1;
j = 2;

while (j < (length(highpass_filt_1)-1) && P_index < R_index)
    if (highpass_filt_1(j-1) > 0 && highpass_filt_1(j+1) <0) && (exp_data_filtered(j)) > .03 && j > R_peaks(P_index)-120,
        output_p(j) = 1; % P-peak
        P_peaks(P_index) = j;
        j = j + 50; % Skip ahead to avoid double counting
        P_index = P_index + 1;
    else
        output_p(j) = 0;
        j = j + 1;
    end
end

figure;
subplot(211);
plot(t1,exp_data,'-.',t1,highpass_filt_3,'-'), grid on
legend('Original Data', 'Scale 3', 2);
subplot(212);
plot(t1,exp_data_filtered,'-.',t1,output_p,'-'), grid on
legend('Data', 'Peaks', 2);

%% 4th scale onward does not show anything meaningful





