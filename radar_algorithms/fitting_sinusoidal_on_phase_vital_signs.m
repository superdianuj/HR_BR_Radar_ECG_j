clear all
close all
clc



%Zhangt, Ting, et al. "Non-contact estimation at 60 GHz for human vital signs monitoring using a robust optimization algorithm." 2016
%IEEE International Symposium on Antennas and Propagation (APSURSI). IEEE, 2016.

%Bakhtiari, Sasan, et al. "A real-time heart rate analysis for a remote millimeter wave IQ sensor."
% IEEE transactions on biomedical engineering 58.6 (2011): 1839-1845.

%Xu, Yu, Qi Li, and Zhenzhou Tang. "Accurate and Contactless Vital Sign Detection in Short Time Window with 24 GHz Doppler Radar."
% Journal of Sensors 2021 (2021).


f1_record=[];
f2_record=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Self-Gathered Experimental Data
for i=1:2

near_pos= table2array(readtable('concentratum_data_50frames_4.csv'));


iChannel=near_pos(:,2);
qChannel=near_pos(:,3);
t=near_pos(:,1);


Fs=1/(t(2)-t(1));
numSecondsBeginning = 1; %Number of seconds to eliminate from beginning of signal
numSecondsEnd = 1;       %Number of seconds to eliminate from end of signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Configuration Details
fileNum = 2;      %2-5

cutoffFreq = 3;          %Highest Frequency to display (Hz)
fPassResp = .2;          %Beginning of passband for respiration rate (Hz)
fStopResp = .5;          %End of passpand for respiration rate (Hz)
fPassHeart = 1;          %Beginning of passband for heart rate (Hz)
fStopHeart = 1.9;        %End of passband for heart rate (Hz)
combWidth = .05;         %width of band to cancel in comb filter
numHarmonics =5;


combinedSignals = iChannel + 1j.*qChannel;
oner=ones(length(iChannel),1);

fun = @(x)sum((abs(iChannel-x(1)).^2+abs(qChannel-x(2)).^2-x(3)*oner.^2).^2);
x0 = [0,0,0];
x = fminsearch(fun,x0);

iChannelp=iChannel-x(1)*oner;
qChannelp=qChannel-x(2)*oner;


theter=atan2(qChannelp,iChannelp);
unwrapped_theter=unwrap(theter);
D_unwrapped_theter=[diff(unwrapped_theter)./diff(t);0];
ampler=qChannelp.^2+iChannelp.^2;

order=4;
framelen=501;
if i==1
    main_sig =  sgolayfilt(theter,order,framelen); % Defining the signal to proceed the algorithm
else
    main_sig =  sgolayfilt(ampler,order,framelen);
end


tries=5;
par_len=7;
P_record=zeros(tries,par_len);
F_record=zeros(tries,1);

for iter=1:tries

    x0=[max(main_sig)*rand
        (fPassResp+fStopHeart)/2*rand
        -100*max(main_sig)/2+100*max(main_sig)*rand
        max(main_sig)*rand
        (fPassHeart+fStopResp)/2*rand
        -100*max(main_sig)/2+100*max(main_sig)*rand
        -max(main_sig)/2+max(main_sig)*rand
        ];

obj= @(x)sum((main_sig-x(1)*sin(x(2)*t+x(3)) ...
        -x(4)*sin(x(5)*t+x(6))-x(7)).^2);

A = [];
b = [];
Aeq = [];
beq = [];

lb = [-max(main_sig),2*pi*fPassHeart,-100*max(main_sig),-max(main_sig),2*pi*fPassResp,-100*max(main_sig),-max(main_sig)];
ub = [max(main_sig),2*pi*fStopHeart, 100*max(main_sig), max(main_sig), 2*pi*fStopResp,100*max(main_sig),max(main_sig)];
options = optimoptions(@fmincon, 'Algorithm' , 'interior-point', 'Display' , 'off');
options.Algorithm = 'sqp' ;
opts.Algorithm = 'interior-point-legacy' ;

[x,fval]=fmincon(obj,x0,A,b,Aeq,beq,lb,ub,[],options);

P_record(iter,:)=x;
F_record(iter)=fval;
end

[M,ind] = min(F_record);
p_estim=P_record(ind,:);



frequency1=p_estim(5)/(4*pi);
frequency2=p_estim(2)/(2*pi);


f1_record=[f1_record;frequency1];
f2_record=[f2_record;frequency2];
end
%% Print out heart and respiration rates

frequency1=min(f1_record);
frequency2=max(f2_record);

endMessage1 = ['Heart Rate is ' num2str(frequency2*60) ...
    ' beats per minute by fitting sinusoidals'];
endMessage2 = ['Respiration Rate is ' num2str(frequency1*60) ...
' breaths per minute by fitting sinusoidals'];
disp(endMessage1);
disp(endMessage2);