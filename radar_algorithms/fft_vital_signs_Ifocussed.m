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
near_pos= table2array(readtable('concentratum_data_50frames.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);



Fs=1/(t(2)-t(1));


%% Configuration Details
fileNum = 2;      %2-5
    if i==1
    inter=2;                 % for breathing rate
    else
        inter=1;                 % for heart rate
    end
order=4;
framelen=1001;
cutoffFreq = 3;          %Highest Frequency to display (Hz)
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


iChannel= sgolayfilt(iChannel,order,framelen);

x_formulized=iChannel;


L = length(x_formulized);   %Length of signals
NFFT = 2^nextpow2(L);   %Length of FFT



%% Take one sided FFT
fftx = fft(x_formulized,NFFT)/L;            
f = Fs/2*linspace(0,1,NFFT/2+1);            %Frequency Range
oneSidedxDFT = 2*abs(fftx(1:NFFT/2+1));


%% Only display frequencies greater than the cutoff frequency
maskCutoff = f>cutoffFreq;
f(maskCutoff) = [];
oneSidedxDFT(maskCutoff) = [];


%% Bandpass filter for respiration rate
respMask = f>fPassResp & f<fStopResp;
xChannelRespDFT = oneSidedxDFT;
xChannelRespDFT(~respMask) = 0;

%% Determine Respiration Rate
[maxxResp , xRespLoc] = max(xChannelRespDFT);
respirationRate = f(xRespLoc);

R_rec=[R_rec, respirationRate];


%% Bandpass filter for heart rate
heartMask = f>fPassHeart & f<fStopHeart;
xChannelHeartDFT = oneSidedxDFT;
xChannelHeartDFT(~heartMask) = 0;

%% Comb filter to eliminate respiration Harmonics
for n = 1:numHarmonics
    combMask = (f < (n*respirationRate + combWidth)) & ...
               (f > (n*respirationRate - combWidth));
    xChannelHeartDFT(combMask) = 0;
end

%% Determine Heart Rate
[maxxHeart , xHeartLoc] = max(xChannelHeartDFT);

heartRate = f(xHeartLoc);
H_rec=[H_rec, heartRate];
end


HR=max(H_rec);
BR=max(R_rec)/2;




%% Print out heart and respiration rates
endMessage1 = ['Heart Rate is ' num2str(HR*60)];
endMessage2 = ['Respiration Rate is ' num2str(BR*60)];
disp(endMessage1);
disp(endMessage2);