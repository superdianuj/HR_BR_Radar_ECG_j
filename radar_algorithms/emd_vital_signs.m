clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data

% Sun, Li, et al. "Remote measurement of human vital signs based on joint-range adaptive EEMD." 
% IEEE Access 8 (2020): 68514-68524.

R_rec=[];
H_rec=[];

near_pos= table2array(readtable('concentratum_data_50frames_4.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);


numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal


%% Configuration Details
fileNum = 2;      %2-5

inter=1;                 % for heart rate


Fs=500*inter;
order=4;
framelen=701;
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

theter=atan2(qChannel,iChannel);
unwrapped_theter=unwrap(theter);

[imf,residual]=emd(unwrapped_theter, 'MaxNumIMF' ,5,'Interpolation' , 'pchip');
rr_peak = 0;
rr_peak1= 0;

max_hr_freq=1.9;
max_br_freq=0.5;

[pks1,locs1]  = findpeaks(imf(:,1),t, 'MinPeakDistance', 1/max_hr_freq, 'MinPeakProminence', 0);
for i = 2:length(locs1)
    rr_peak = rr_peak + locs1(i) - locs1(i-1);
end

HR = ((length(locs1)-1)/rr_peak)*60;


[pks,locs]  = findpeaks(imf(:,3),t, 'MinPeakDistance', 1/max_br_freq, 'MinPeakProminence', 0);
for i = 2:length(locs)
    rr_peak1 = rr_peak1 + locs(i) - locs(i-1);
end

RR = ((length(locs)-1)/rr_peak)*60/2;

%% Print out heart and respiration rates
endMessage1 = ['Heart Rate is ' num2str(HR) ...
    ' beats per minute'];
endMessage2 = ['Respiration Rate is ' num2str(RR) ...
' breaths per minute'];
disp(endMessage1);
disp(endMessage2);

