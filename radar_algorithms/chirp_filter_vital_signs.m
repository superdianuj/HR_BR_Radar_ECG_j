clear all
close all
clc

% Nejadgholi, I., et al. "Estimation of breathing rate with confidence interval using single-channel CW radar."
% Journal of Healthcare Engineering 2019 (2019).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data
near_pos= table2array(readtable('concentratum_data_50frames.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);
order=4;
framelen=931;
iChannel= sgolayfilt(iChannel,order,framelen);
qChannel= sgolayfilt(qChannel,order,framelen);

Fs=1/(t(2)-t(1));
numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(t,iChannel)
hold on
plot(t,qChannel)
grid on
xlabel('time(s)')
ylabel('Volatage (mV)')


%% Configuration Details

cutoffFreq = 3;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .5;          %End of passpand for respiration rate (Hz)
fPassHeart = 1;          %Beginning of passband for heart rate (Hz)
fStopHeart = 1.8;        %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics =5;         %number of harmonics to cancel in comb filter


combinedSignals = iChannel + 1j.*qChannel;
oner=ones(length(iChannel),1);


fun = @(x)sum((abs(iChannel-x(1)).^2+abs(qChannel-x(2)).^2-x(3)*oner.^2).^2);
x0 = [0,0,0];
x = fminsearch(fun,x0);


iChannelp=iChannel-x(1)*oner;
qChannelp=qChannel-x(2)*oner;


theter=atan2(qChannelp,iChannelp);
unwrapped_theter=unwrap(theter);

ampler=qChannelp.^2+iChannelp.^2;
Signal_K = detrend((unwrapped_theter-mean(unwrapped_theter)));  %deterend the signal

[b,a] = butter(5,10/Fs,'low'); % 5Hz lowpass
Signal_K = filter(b,a,(Signal_K));

Signal_K_t = Signal_K;  % keep a copy of bandpass filetered signal


x_br_fft= chirp_based_estimator( Signal_K,Fs ,fPassResp,fStopResp);  % Estiamtion of breathing from radar using Chirp transform
x_br_fft = 60*x_br_fft/2;  


x_hr_fft= chirp_based_estimator( Signal_K,Fs,fPassHeart,fStopHeart);  % Estiamtion of breathing from radar using Chirp transform
x_hr_fft = 60*x_hr_fft;  

endMessage1 = ['Breathing Rate is ' num2str(x_br_fft) ...
    ' beats per minute'];
disp(endMessage1);
endMessage1 = ['Heart Rate is ' num2str(x_hr_fft) ...
    ' beats per minute'];
disp(endMessage1);