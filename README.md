Sup guys.

This is the readme file for the repository. We're going to be attemtping to write up an algo to analyze this ECG signal for AF. I've added a script here called
fft.m and it takes a dataset on a text document and runs an FFT on it. This is from the script I used to test a microphone, so there are some extra parts that 
we can remove (the part where I recreate the sound on a .wav file). 

So in order to recreate this algorithm, we'll need to read the papers that outline how to do this.
I'll find the papers that do so and post them here:

Wavelet Transform:
-Instead of decomposing into sin/cos you instead break signal into wavelets, all based on a mother wavelet
-If mother wavelet is derivative of smoothing function, then WT will be proportional to derivative of signal
-Instead of using integrals, implement a filter bank which replicates a Discrete WT