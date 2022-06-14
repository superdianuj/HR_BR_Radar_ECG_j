clear all
close all
clc

%Seflek, Ibrahim, Yunus Emre Acar, and Ercan Yaldiz. "Small motion detection and non-contact vital signs monitoring with continuous wave Doppler radars." 
%Elektronika ir Elektrotechnika 26.3 (2020): 54-60.

%Yang, Zi-Kai, et al. "Accurate Doppler radar-based heart rate measurement using matched filter."
% IEICE Electronics Express 17.8 (2020): 20200062-20200062.

%Srihari, Pathipati, and G. S. Vandana. "Experimental Study of 24GHz Sense2Gol Pulse Radar Sensor for Human Vital Sign Measurement." 2021 
%IEEE International Conference on Electronics, Computing and Communication Technologies (CONECCT). IEEE, 2021.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data

R_rec=[];
H_rec=[];
for i=1:2
near_pos= table2array(readtable('concentratum_data_50frames_1.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);


numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal




%% Configuration Details
fileNum = 2;      %2-5
    if i==1
    inter=2;                 % for breathing rate
    else
        inter=1;                 % for heart rate
    end

Fs=500*inter;
order=4;
framelen=501;
cutoffFreq = 5;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .5;          %End of passpand for respiration rate (Hz)
fPassHeart = 1;          %Beginning of passband for heart rate (Hz)
fStopHeart = 1.8;        %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics =5;  



t_new=linspace(0,max(t),length(t)*inter);


IC=spline(t,iChannel,t_new);
QC=spline(t,qChannel,t_new);

iChannel=IC';
qChannel=QC';
t=t_new';



oner=ones(length(iChannel),1);
fun = @(x)sum((abs(iChannel-x(1)).^2+abs(qChannel-x(2)).^2-x(3)*oner.^2).^2);
x0 = [0,0,0];
x = fminsearch(fun,x0);


iChannel=iChannel-x(1)*oner;
qChannel=qChannel-x(2)*oner;


iChannel= sgolayfilt(iChannel,order,framelen);
qChannel= sgolayfilt(qChannel,order,framelen);


combinedSignals = iChannel + 1j.*qChannel;

% theter=atan2(qChannel,iChannel);
% unwrapped_theter=unwrap(theter);
% 
% D_unwrapped_theter=[diff(unwrapped_theter);0];
% 
% combinedSignals=unwrapped_theter;


L = length(iChannel);   %Length of signals
NFFT = 2^nextpow2(L);   %Length of FFT



%% Take one sided FFT
fftI = fft(iChannel,NFFT)/L;                %FFT of I channel
fftQ = fft(qChannel,NFFT)/L;                %FFT of Q channel
fftCombined = fft(combinedSignals,NFFT)/L;  %FFT of Q channel
f = Fs/2*linspace(0,1,NFFT/2+1);            %Frequency Range
oneSidedIDFT = 2*abs(fftI(1:NFFT/2+1));
oneSidedQDFT = 2*abs(fftQ(1:NFFT/2+1));
oneSidedCombinedDFT = 2*abs(fftCombined(1:NFFT/2+1));

%% Only display frequencies greater than the cutoff frequency
maskCutoff = f>cutoffFreq;
f(maskCutoff) = [];
oneSidedIDFT(maskCutoff) = [];
oneSidedQDFT(maskCutoff) = [];
oneSidedCombinedDFT(maskCutoff) = [];

%% Bandpass filter for respiration rate
respMask = f>fPassResp & f<fStopResp;
iChannelRespDFT = oneSidedIDFT;
qChannelRespDFT = oneSidedQDFT;
combinedRespDFT = oneSidedCombinedDFT;
iChannelRespDFT(~respMask) = 0;
qChannelRespDFT(~respMask) = 0;
combinedRespDFT(~respMask) = 0;

%% Determine Respiration Rate
[maxIResp , iRespLoc] = max(iChannelRespDFT);
[maxQResp , qRespLoc] = max(qChannelRespDFT);
[maxCombinedResp , combinedRespLoc] = max(combinedRespDFT);

respirationRate = f(combinedRespLoc);
respChoice = 'Combined Channel';

if(maxIResp > maxQResp && maxIResp > maxCombinedResp)
    respirationRate = f(iRespLoc);
    respChoice = 'I channel';
end
if(maxQResp > maxIResp && maxQResp > maxCombinedResp)
    respirationRate = f(qRespLoc);
    respChoice = 'Q channel';
end
if(maxCombinedResp > maxIResp && maxCombinedResp> maxQResp)
    respirationRate = f(combinedRespLoc);
    respChoice = 'Combined channel';
end

R_rec=[R_rec, respirationRate];
%% Bandpass filter for heart rate
heartMask = f>fPassHeart & f<fStopHeart;
iChannelHeartDFT = oneSidedIDFT;
qChannelHeartDFT = oneSidedQDFT;
combinedHeartDFT = oneSidedQDFT;
iChannelHeartDFT(~heartMask) = 0;
qChannelHeartDFT(~heartMask) = 0;
combinedHeartDFT(~heartMask) = 0;

%% Comb filter to eliminate respiration Harmonics
for n = 1:numHarmonics
    combMask = (f < (n*respirationRate + combWidth)) & ...
               (f > (n*respirationRate - combWidth));
    iChannelHeartDFT(combMask) = 0;
    qChannelHeartDFT(combMask) = 0;
    combinedHeartDFT(combMask) = 0;
end

%% Determine Heart Rate
[maxIHeart , iHeartLoc] = max(iChannelHeartDFT);
[maxQHeart , qHeartLoc] = max(qChannelHeartDFT);
[maxCombinedHeart , combinedHeartLoc] = max(combinedHeartDFT);

heartRate = f(combinedHeartLoc);
heartChoice = 'Combined Channel';

if(maxIHeart > maxQHeart && maxIHeart > maxCombinedHeart)
    heartRate = f(iHeartLoc);
    heartChoice = 'I channel';
end
if(maxQHeart > maxIHeart && maxQHeart > maxCombinedHeart)
    heartRate = f(qHeartLoc);
    heartChoice = 'Q channel';
end
if(maxCombinedHeart > maxIHeart && maxCombinedHeart > maxQHeart)
    heartRate = f(combinedHeartLoc);
    heartChoice = 'Combined channels';
end
H_rec=[H_rec, heartRate];
end


HR=max(H_rec);
BR=max(R_rec)/2;
%% Plot I and Q
figure
subplot(3,1,1)
plot(t,iChannel)
xlabel('Time (s)')
ylabel('|i(t)|')
title('I Channel in Time Domain')
subplot(3,1,2)
plot(t,qChannel)
xlabel('Time (s)')
ylabel('|q(t)|')
title('Q Channel in Time Domain')
subplot(3,1,3)
plot(t,abs(combinedSignals))
xlabel('Time (s)')
ylabel('|c(t)|')
title('Combined Signals in Time Domain')

figure
subplot(3,1,1)
plot(f,oneSidedIDFT) 
title('Single-Sided Amplitude Spectrum of I channel FFT')
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
subplot(3,1,2)
plot(f,oneSidedQDFT) 
title('Single-Sided Amplitude Spectrum of Q channel FFT')
xlabel('Frequency (Hz)')
ylabel('|Q(f)|')
subplot(3,1,3)
plot(f,oneSidedCombinedDFT) 
title('Single-Sided Amplitude Spectrum of Combined channels FFT')
xlabel('Frequency (Hz)')
ylabel('|C(f)|')

figure
subplot(3,1,1)
plot(f,iChannelRespDFT) 
title('I Channel Bandpass for Respiration Rate')
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
subplot(3,1,2)
plot(f,qChannelRespDFT) 
title('Q Channel Bandpass for Respiration Rate')
xlabel('Frequency (Hz)')
ylabel('|Q(f)|')
subplot(3,1,3)
plot(f,combinedRespDFT) 
title('Combined Channels Bandpass for Respiration Rate')
xlabel('Frequency (Hz)')
ylabel('|C(f)|')

figure
subplot(3,1,1)
plot(f,iChannelHeartDFT) 
title('I Channel Bandpass for Heart Rate')
xlabel('Frequency (Hz)')
ylabel('|I(f)|')
subplot(3,1,2)
plot(f,qChannelHeartDFT) 
title('Q Channel Bandpass for Heart Rate')
xlabel('Frequency (Hz)')
ylabel('|Q(f)|')
subplot(3,1,3)
plot(f,combinedHeartDFT) 
title('Combined Channels Bandpass for Heart Rate')
xlabel('Frequency (Hz)')
ylabel('|C(f)|')

%% Print out heart and respiration rates
endMessage1 = ['Heart Rate is ' num2str(HR*60) ...
    ' beats per minute using the ' heartChoice];
endMessage2 = ['Respiration Rate is ' num2str(BR*60) ...
    ' breaths per minute using the ' respChoice];
disp(endMessage1);
disp(endMessage2);