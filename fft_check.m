function [ ] = fft_check( Fs, data, title_graph )

% Define the axis from which calculate the fft
axis = 1;
% Define the stabilization period (percentage)
Stab = 0.1;
% Calculate the period
T = 1/Fs;
% Get the size of the dataset
size = length(data(:,axis));
% Calculate a stabilization period
stabilization = round(Stab * size);
% Remove the data from the beginning and end of the dataset
data = data(stabilization+1:(size-stabilization),:);
% Calculate the new size of the dataset 
size = length(data(:,axis));
% Creates a vector for the time domain
t = (0:size-1)*T;
% Calculate the next power of two of the dataset (speed up FFT)
NFFT = 2^nextpow2(size);
% Calculate the fft for the dataset, for NFFT points
fftData =(abs(fft(data,NFFT)));
%shift to zero frequency
fftshift(fftData);
% Calculate the frequency vector
freq = Fs*(0:NFFT/2-1)/NFFT;
% Plot the frequency spectrum
figure;
subplot(2,1,1);
plot(freq,fftData(1:NFFT/2, axis));
title(title_graph);
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
% Plot the original data
subplot(2,1,2);
plot(t,data(:, axis)) ;
title('Original Data');
xlabel('Time (s)');
ylabel('Amplitude');
end


