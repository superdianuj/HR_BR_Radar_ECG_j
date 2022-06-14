clear all
close all
clc


%Zheng, Tianyue, et al. "V2ifi: In-vehicle vital sign monitoring via compact rf sensing."
%Proceedings of the ACM on Interactive, Mobile, Wearable and Ubiquitous Technologies 4.2 (2020): 1-27.

%Zheng, Tianyue, et al. "V2ifi: In-vehicle vital sign monitoring via compact rf sensing."
%Proceedings of the ACM on Interactive, Mobile, Wearable and Ubiquitous Technologies 4.2 (2020): 1-27.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data
near_pos= table2array(readtable('concentratum_data_50frames_4.csv'));


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


%% Configuration Details
fileNum = 2;      %2-5

cutoffFreq = 5;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .9;          %End of passpand for respiration rate (Hz)
fPassHeart = 2;          %Beginning of passband for heart rate (Hz)
fStopHeart = 3;          %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics = 5;        %number of harmonics to cancel in comb filter


combinedSignals = iChannel + 1j.*qChannel;
oner=ones(length(iChannel),1);


fun = @(x)sum((abs(iChannel-x(1)).^2+abs(qChannel-x(2)).^2-x(3)*oner.^2).^2);
x0 = [0,0,0];
x = fminsearch(fun,x0);


iChannelp=iChannel-x(1)*oner;
qChannelp=qChannel-x(2)*oner;


ampler=iChannelp.^2+qChannelp.^2;

ampler_squarred=2*ampler.^2;

interpolation_factor=1;
decimation_factor=2;
filty=designMultirateFIR(interpolation_factor,decimation_factor);
firDecim = dsp.FIRDecimator(decimation_factor,filty);
siggy = firDecim(ampler(1:end));

desired_sam_factor=2;
siggy_down=downsample(siggy,desired_sam_factor);

passband_freq=1;
new_sampling_frequency=length(siggy_down)/(t(end-1)-t(1));
result_sig=lowpass(siggy_down,passband_freq,new_sampling_frequency);
envelope=abs(result_sig).^0.5;


timer=linspace(t(1),t(end-1),length(envelope));
env = [];
rr_peak = 0;
rr_peak1=0;

if mod(length(envelope),2) == 0 %Ensuring l is always even as FFT requires an even length
    l = length(envelope);
else
    l = length(envelope)-1;
end
    
env=envelope;
t = timer';

imf  = vmd(env,'NumIMFs',4);

max_hr_freq=130/60;
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