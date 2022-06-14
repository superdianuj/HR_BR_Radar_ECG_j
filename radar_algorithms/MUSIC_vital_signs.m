clear all
close all
clc

% Nejadgholi, Isar, Sreeraman Rajan, and Miodrag Bolic. "Time-frequency based contactless estimation of vital signs of human while walking using PMCW radar." 2016 IEEE 18th International Conference on e-Health Networking, 
% Applications and Services (Healthcom). IEEE, 2016.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data
near_pos= table2array(readtable('concentratum_data_50frames_2.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);
inter=1;
t_new=linspace(0,max(t),length(t)*inter);


IC=spline(t,iChannel,t_new);
QC=spline(t,qChannel,t_new);

iChannel=IC';
qChannel=QC';
t=t_new';

Fs=1/(t(2)-t(1));

numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Configuration Details

cutoffFreq = 3;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .5;          %End of passpand for respiration rate (Hz)
fPassHeart = 1;          %Beginning of passband for heart rate (Hz)
fStopHeart = 1.8;        %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics =5;         %number of harmonics to cancel in comb filter




ampler=qChannel.^2+iChannel.^2;


% normalizing
Signal_thisZone=ampler;

Signal = (Signal_thisZone-mean(Signal_thisZone))/(max(Signal_thisZone) - min(Signal_thisZone));

L = size(Signal,1);
t =(0:(L-1))/Fs;
TargetFmax_matrix = [];  %Matrix of Observations to be filled for this 20 seconds


St = Signal;
L_St = length(St);
t_St = (0:(L_St-1))/Fs;

% calculating TFR of the detected signal
%
WT =[];
freq = [];
wopt = [];
[WT,freq,wopt]=stft(St,Fs,'Window' ,kaiser(256,5), 'OverlapLength' ,220, 'FFTLength' ,512);

%%% Find instantanious frequencies related to maximum amplitude
Fmax = [];
for nn = 1:length(wopt)
    [a,b(nn)]= max(abs(WT(:,nn)));
    Fmax(nn) = freq(b(nn));
end
Fmax = Fmax-mean(Fmax);  % Remove DC

%add the observed instantanious frequency to the Observation matrix
TargetFmax_matrix= [TargetFmax_matrix ;Fmax];


% calculate MUSIC spectrum
[S,f]= pmusic(TargetFmax_matrix,4,2^16,Fs);

% Finding Breathing rate

%% Bandpass filter for respiration rate
respMask = f>fPassResp & f<fStopResp;
xChannelRespDFT = S;
xChannelRespDFT(~respMask) = 0;

%% Determine Respiration Rate
[maxxResp , xRespLoc] = max(xChannelRespDFT);
Breathing_rate = f(xRespLoc)*60;
respirationRate=Breathing_rate/(60*2);


%% Bandpass filter for heart rate
heartMask = f>fPassHeart & f<fStopHeart;
xChannelHeartDFT = S;
xChannelHeartDFT(~heartMask) = 0;

%% Comb filter to eliminate respiration Harmonics
for n = 1:numHarmonics
    combMask = (f < (n*respirationRate + combWidth)) & ...
               (f > (n*respirationRate - combWidth));
    xChannelHeartDFT(combMask) = 0;
end

%% Determine Heart Rate
[maxxHeart , xHeartLoc] = max(xChannelHeartDFT);

Heart_rate = f(xHeartLoc)*60;


%% Print out heart and respiration rates
endMessage1 = ['Heart Rate is ' num2str(Heart_rate) ...
    ' beats per minute'];
endMessage2 = ['Respiration Rate is ' num2str(Breathing_rate) ...
' breaths per minute'];
disp(endMessage1);
disp(endMessage2);