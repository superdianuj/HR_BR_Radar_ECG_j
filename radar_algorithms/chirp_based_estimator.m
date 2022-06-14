function [x_fr] = chirp_based_estimator( Signal_K,fs,low_freq,high_freq)
Data_RS = Signal_K;

[b,a] = butter(4,0.1/fs,'high'); % 0.175Hz highpass
Data_RS = filter(b,a,detrend(Data_RS));
[b,a] = butter(5,2/fs,'low'); % 4Hz lowpass
Data_RS = filter(b,a,(Data_RS));

Data_RS = Data_RS(:).*hamming(length(Data_RS));
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Chirp Transform %%%%%%%%%

%% Search for Breathing Rate

m=2048; % number of sample points to calculate for chirp transform
f1 = 0; % lower frequency bound of chirp transform
f2 = 4; % upper frequency bound of chirp transform
w = exp(-j*2*pi*(f2-f1)/(m*fs)); % arc in unit circle of z domain is defined by w and a
a = exp(j*2*pi*f1/fs);
z = czt(Data_RS,m,w,a); % chirp transform
fn = (0:m-1)'/m; % normalized frequency vector
fy = fs*fn; % un-normalized frequency vector
fz = (f2-f1)*fn + f1; % adding back f1 (in case f1 is not zero)
freqIndex = find(fz>low_freq & fz<high_freq );
[Ma I]=max(abs(z(freqIndex)));
x_fr= fz(I+freqIndex(1));  
end