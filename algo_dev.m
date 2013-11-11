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

%%% ECG Data from MIT Database %%%
[t1,exp_data] = rdsamp('mitdb/102', 1,3800,200); % 3600 points which corresponds to 10 s @ 360 sampling rate

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
b1 = [1/8, 3/8, 0, 3/8];

lowpass_filt_1 = filter(b1,a1,exp_data_filtered);

%Plot to see how it looks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(t,exp_data,'-.',t,lowpass_filt,'-'), grid on
%legend('Original Data', 'Low-passed Data', 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Design FIR filter as 2 * (d[n+1] - d[n])

a2 = 1;
b2 = [0, 2, -2, 0];

highpass_filt = filter(b2,a2,exp_data);
highpass_filt_1 = filter(b2,a2,exp_data_filtered);

output = zeros(size(highpass_filt_1));
lengths = zeros(size(highpass_filt_1));
RRmean = zeros(size(highpass_filt_1));

R_peaks = zeros(300:1); % Assume max heart rate is 300 bpm and we're going to be making classification every minute
R_index = 1;
P_counter = 1;
j = 2;
interval_sum = 0;

% Pick out Important Maxima (Classification)%
while (j < (length(highpass_filt_1)-1))
    if (highpass_filt_1(j-1) > 0 && highpass_filt_1(j+1) <0) && (exp_data_filtered(j)) > .5,
        output(j) = 1; % R-peak
        R_peaks(R_index) = j;
        if(R_index == 1)
            continue;
        else
            lengths(R_index-1) = R_peaks(R_index) - R_peaks(R_index - 1);
            interval_sum = interval_sum + lengths(R_index - 1);
        end
        j = j + 50; % Skip ahead to avoid double counting
    else
        output(j) = 0;
        j = j + 1;
    end
end

%First Scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% t2 = 1:length(first_scale);
% subplot(211);
% plot(t1,exp_data_filtered,'-.',t1,highpass_filt_1,'-'), grid on
% legend('Data', 'Scale-1 Data', 2);
% subplot(212);
plot(t1,exp_data_filtered,'-.',t1,output,'-'), grid on
legend('Data', 'Peaks', 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Interpolate signal by 2 before being fed into 2nd scale filters%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Low pass Butterworth filter of order 10 at Nyquist Frequency
[b_inter,a_inter] = butter(10,.99,'low'); 

% First thing to do is upsample by a factor of 2
y = upsample(lowpass_filt_1,2);

% Then run it through a low-pass filter at half the sampling freq
lowpass_filt_inter = filtfilt(b_inter,a_inter, y); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %2nd Scale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% highpass_filt_2 = filter(b2,a2,lowpass_filt_inter);
% lowpass_filt_2 = filter(b1,a1, lowpass_filt_inter);
% 
% % To plot highpass_filt_2 output, I downsample by 2
% highpass_filt_2 = downsample(highpass_filt_2, 2);
% 
% figure;
% plot(t1,exp_data_filtered,'-.',t1,highpass_filt_2,'-'), grid on
% legend('Original Data', 'Scale 2', 2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %3rd Scale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% highpass_filt_3 = filter(b2,a2,lowpass_filt_2);
% lowpass_filt_3 = filter(b1,a1, lowpass_filt_2);
% figure;
% plot(t1,exp_data,'-.',t1,highpass_filt_3,'-'), grid on
% legend('Original Data', 'Scale 3', 2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %4th Scale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% highpass_filt_4 = filter(b2,a2,lowpass_filt_3);
% lowpass_filt_4 = filter(b1,a1, lowpass_filt_3);
% figure;
% plot(t1,exp_data,'-.',t1,highpass_filt_4,'-'), grid on
% legend('Original Data', 'Scale 4', 2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %5th Scale
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% highpass_filt_5 = filter(b2,a2,lowpass_filt_4);
% lowpass_filt_5 = filter(b1,a1, lowpass_filt_4);
% figure;
% plot(t1,exp_data,'-.',t1,highpass_filt_5,'-'), grid on
% legend('Original Data', 'Scale 5', 2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Find R peaks%%%





