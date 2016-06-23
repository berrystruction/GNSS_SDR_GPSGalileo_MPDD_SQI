clc
clear all
close all
addpath('include','geoFunctions','E:\DataSet\MPMP\newModelMP\5MPtimeWin\5MPTime_4sat')

% Pseudorate Values available
% 250, 301, 375, 400, 500, 1000, 1044 ms (pseudorange output) 5 sat
% 103, 301, 500, 1000, 2000 ms (pseudorange output) 4 sat
pseudorate=500;
%load(['trackingResults_data_CN0_100_NoMP_2bit_60s_2_pseudo_' num2str(pseudorate) 'ms' ]);
%load trackingResults_data_pseudo500ms_16368_41304_noMP_NONOISE_NOSKIP;
load 'trackingResults_data_16368_41304_noMP_NONOISE_pseudo_500ms';
settings1=settings;
trackResults1=trackResults;
navSolutions1 = postNavigation(trackResults1, settings1);

clearvars -except navSolutions1 settings1 trackResults1 pseudorate;

%load(['trackingResults_data_ant_1_noise_AllMP_2bit_60s_5diffWin_pseudo_' num2str(pseudorate) 'ms' ]);
load(['trackingResults_data_CN0_100_NoMP_2bit_60s_2_pseudo_' num2str(pseudorate) 'ms' ]);
%load 'trackingResults_data_16368_41304_CN0_100_NoMP_2bit_60s_pseudo_500ms_elevmask.mat'
%load(['trackingResults_data_CN0_100_NoMP_2bit_60s_pseudo_' num2str(pseudorate) 'ms' ]);
%load 'trackingResults_data_16368_41304_noMP_NONOISE_pseudo_500ms_SKIP3s';
%trackingResults_data_pseudo1000ms_16368_41304_noMP_NONOISE_NOSKIP;
%trackingResults_data_pseudo1000ms_noise_AllMP_2bit_60s_differentWin
settings2=settings;
trackResults2=trackResults;

navSolutions2 = postNavigation(trackResults2, settings2);
%%%%%%%%%%%%%%%%%%%%%
plotNavigation(navSolutions1, settings1);
plotNavigation(navSolutions2, settings2);

discreteTime=[0:length(navSolutions2.channel.rawP(end,:))-1];
xtime=(discreteTime+...
    floor(settings.msToProcess/settings.navSolPeriod-length(navSolutions2.channel.rawP(end,:))))*settings.navSolPeriod/1e3;

load ('E:\DataSet\MPMP\newModelMP\5MPtimeWin\params','MPS')
%%
timeMP=[];
for indexperiod=1:length(MPS)
    timeMP=[timeMP MPS(indexperiod).TimePeriod];
end
timeMP=timeMP/1e3;

profileMP1=zeros(size(discreteTime));
index1=1;
index2=2;
for indexperiod= 1:length(discreteTime)
    if xtime(indexperiod)>=timeMP(index1) && xtime(indexperiod)<=timeMP(index2)
        profileMP1(indexperiod)=1;
    end
    if xtime(indexperiod)>timeMP(index2) && index2<length(timeMP)
        index1=index1+2;
        index2=index2+2;
    end
end
%%
millisec=1;
timelimit=floor(xtime(end)/(millisec*1e-3));
profileMP2=zeros(1,timelimit);
index1=1;
index2=2;
timeMP=timeMP*1e3;
for indexperiod= 1:timelimit
    if indexperiod*millisec>=timeMP(index1) && (indexperiod*millisec)<=timeMP(index2)
        profileMP2(indexperiod)=1;
    end
    if (indexperiod*millisec)>timeMP(index2) && index2<length(timeMP)
        index1=index1+2;
        index2=index2+2;
    end
end
ttime=[1:timelimit]*millisec*1e-3;

%%
for index=1:size(navSolutions2.channel.PRN,1)
    %[neworder1 I1]=sort(navSolutions1.channel.PRN(:,1));
    %     [neworder2 I2]=sort(navSolutions2.channel.PRN(:,1));
    %     plot(navSolutions2.channel.rawP(I2(index),1:end),'.-')
    
    figure
    subplot(221)
    [neworder1 I1]=sort(navSolutions1.channel.PRN(:,1));
    [neworder2 I2]=sort(navSolutions2.channel.PRN(:,1));
    %I2=[5 2 4 3 1];%9 15 18 12 14];
    
    plot(xtime,navSolutions1.channel.rawP(I1(index),:)-navSolutions2.channel.rawP(I2(index),1:end),'.-',...
    xtime,profileMP1*max(abs(navSolutions1.channel.rawP(I1(index),:)-navSolutions2.channel.rawP(I2(index),1:end))),'ro')
    grid on
    title(['RAW (\rho-\rho_{MP}). PRN:' num2str(navSolutions2.channel.PRN(I2(index),1))])
    xlabel('Time [s]')
    
    
    subplot(222)
    plot(xtime,navSolutions2.channel.rawP(I2(index),:),'.')
    grid on
    title('RAW \rho')
    xlabel('Time [s]')
    
    subplot(223)
    plot(ttime,profileMP2,'r.-')
    grid on
    title('MP profile')
    xlabel('Time [s]')
    
    subplot(224)
    plot(xtime,profileMP1,'r.-')
    grid on
    title(['MP profile (pseudorange rate output: ' num2str(settings.navSolPeriod) ' ms)'])
    xlabel('Time [s]')
    
    
end

figure
plot(xtime,navSolutions1.E-navSolutions2.E,'*-',xtime,navSolutions1.N-navSolutions2.N,'o-',xtime,navSolutions1.U-navSolutions2.U,'^-')
title('Difference in navigation solution')
legend('East','North','Up')
grid on

timeMP=timeMP/1e3;
fprintf('Multipath istants in seconds:\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n',...
    timeMP(1),timeMP(2),timeMP(2)-timeMP(1),...
    timeMP(3),timeMP(4),timeMP(4)-timeMP(3),...
    timeMP(5),timeMP(6),timeMP(6)-timeMP(5),...
    timeMP(7),timeMP(8),timeMP(8)-timeMP(7),...
    timeMP(9),timeMP(10),timeMP(10)-timeMP(9));
%h = msgbox({'multipath istants in seconds:' num2str(timeMP/1e3) },'MP windows');