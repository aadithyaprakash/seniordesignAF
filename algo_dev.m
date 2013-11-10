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

experiment_number = 0;
experiment = ['ecg_test_', num2str(experiment_number), '.txt'];
exp_data = load(experiment);


% Design FIR Filter as (1/8) * (d[n+2]+ 3d[n+1] + 3d[n-1]) where d is dirac
% impulse (NON-CAUSAL so first few ouputs will be garbage)

a1 = 1;
b1 = [1/8, 3/8, 0, 3/8];

lowpass_filt = filter(b1,a1,exp_data);

t = 1:length(exp_data);
plot(t,exp_data,'-.',t,lowpass_filt,'-'), grid on
legend('Original Data', 'Low-passed Data', 2);


% Design FIR filter as 2 * (d[n+1] - d[n])

a2 = 1;
b2 = [0, 2, -2, 0];

highpass_filt = filter(b2,a2,exp_data);
figure;
plot(t,exp_data,'-.',t,highpass_filt,'-'), grid on
legend('Original Data', 'High-passed Data', 2);