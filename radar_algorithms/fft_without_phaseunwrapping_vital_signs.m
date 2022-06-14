clear all
close all
clc

%Srihari, Pathipati, and G. S. Vandana. "Experimental Study of 24GHz Sense2Gol Pulse Radar Sensor for Human Vital Sign Measurement." 2021 
%IEEE International Conference on Electronics, Computing and Communication Technologies (CONECCT). IEEE, 2021.

%Seflek, Ibrahim, Yunus Emre Acar, and Ercan Yaldiz. "Small motion detection and non-contact vital signs monitoring with continuous wave Doppler radars." 
%Elektronika ir Elektrotechnika 26.3 (2020): 54-60.

%Yang, Zi-Kai, et al. "Accurate Doppler radar-based heart rate measurement using matched filter."
% IEICE Electronics Express 17.8 (2020): 20200062-20200062.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data

% Sun, Li, et al. "Remote measurement of human vital signs based on joint-range adaptive EEMD." 
% IEEE Access 8 (2020): 68514-68524.

R_rec=[];
H_rec=[];

near_pos= table2array(readtable('concentratum_data_50frames_2.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);


numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal


%% Configuration Details

inter=1;                


t_new=linspace(0,max(t),length(t)*inter);


IC=spline(t,iChannel,t_new);
QC=spline(t,qChannel,t_new);

iChannel=IC';
qChannel=QC';
t=t_new';


oner=ones(length(iChannel),1);
fun = @(x)sum((iChannel-x(1)/x(2)).^2+(((qChannel-x(3))/(x(4)*cos(x(5)))-((iChannel-x(1))*sin(x(5))/(x(2)*cos(x(5))))).^2-oner.^2).^2);
x0 = [0,0,0,0,0];
x = fminsearch(fun,x0);

dci=x(1);
ai=x(2);
dcq=x(3);
aq=x(4);
phi=x(5);

iChannel=iChannel-dci;
qChannel=(qChannel-dcq-(ai/aq)*iChannel*sin(phi))/cos(phi);

Ic=iChannel(1:end-1);
DQc=diff(qChannel);
Qc=qChannel(1:end-1);
DIc=diff(iChannel);

x_formulized=zeros(length(Ic),1);
for i=1:length(x_formulized)
    ic=Ic(1:i);
    dqc=DQc(1:i);
    qc=Qc(1:i);
    dic=DIc(1:i);
    x_formulized(i)=sum((ic.*dqc-qc.*dic)./((aq/ai)*ic.^2+(ai/aq)*qc.^2));
end


% order=2;
% framelen=101;
% x_formulized = sgolayfilt(x_formulized,order,framelen); % Defining the signal to proceed the algorithm


t=t(1:end-1);

Fs=1/(t(2)-t(1));

%% Configuration Details



cutoffFreq = 3;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .5;          %End of passpand for respiration rate (Hz)
fPassHeart = 1;          %Beginning of passband for heart rate (Hz)
fStopHeart = 1.8;        %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics =5;  



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

%-----------------

HR=heartRate;
BR=respirationRate/2;




%% Print out heart and respiration rates
endMessage1 = ['Heart Rate is ' num2str(HR*60)];
endMessage2 = ['Respiration Rate is ' num2str(BR*60)];
disp(endMessage1);
disp(endMessage2);


