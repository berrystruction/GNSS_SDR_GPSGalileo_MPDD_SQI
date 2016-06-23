clc
clearvars
close all
%addpath('include','geoFunctions','E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\NoNoise\');%nearFarMP\trackresults_widefilter\')
%  'E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\60dBHz\trackresults\T_int_1ms')%'E:\DataSet\MPMP\newModelMP\5MPtimeWin\5MPTime_4sat')%%,'E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph')
% Comparison between data coming from analysis of dataset with different
% pseudorange output rate with/without multipath

% Pseudorate Values available
% 250, 301, 375, 400, 500, 1000, 1044 ms (pseudorange output) 5 sat
%pseudorate=[250 301 375 400 500 1000 1044];
% 103, 301, 500, 1000, 2000 ms (pseudorange output) 4 sat
%pseudorate=[103 301 500 1000 2000];
pseudorate=[500 1000 2000];

Ndataset=length(pseudorate);
startPreamblems=zeros(size(pseudorate));

settingsV0=[];
trackResultsV0=[];
navSolutionsV0=[];
settingsV1=[];
trackResultsV1=[];
navSolutionsV1=[];
settingsV2=[];
trackResultsV2=[];
navSolutionsV2=[];
settingsV3=[];
trackResultsV3=[];
navSolutionsV3=[];
settingsV4=[];
trackResultsV4=[];
navSolutionsV4=[];
settingsV5=[];
trackResultsV5=[];
navSolutionsV5=[];
% Loading data
disable=1;
for datasetindex=1:Ndataset
    
    disp(['loading dataset data #' num2str(datasetindex)]);
    
    %switch datasetindex
        %case 1
        %    load('/Users/mattiaberardo/Dataset Thuan/mytrackingResults9404_trackingResults_data_2_pseudo_500ms_Bl20.mat')
        %case 2
        %    load('/Users/mattiaberardo/Dataset Thuan/mytrackingResults9404_trackingResults_data_2_pseudo_00ms_Bl20.mat')    
        %case 3
        %    load('/Users/mattiaberardo/Dataset Thuan/mytrackingResults9404_trackingResults_data_2_pseudo_500ms_Bl20.mat')         
    %end
    
    %load(['/Volumes/Elements/BackUp_Datasets/DataSets/MPMP/newModelMP/5MPtimeWin/differentEph/2bits/NoNoise/trackingResults_maggio2016/trackingResults_data_0_pseudo_' num2str(pseudorate(datasetindex)) 'ms']); % ideale classico')         
    load(['/Volumes/Elements/BackUp_Datasets/DataSets/MPMP/newModelMP/5MPtimeWin/differentEph/2bits/NoNoise/trackresults/trackingResults_data_0_pseudo_' num2str(pseudorate(datasetindex)) 'ms']); % ideale classico

    %load(['E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\NoNoise\trackresults\trackingResults_data_0_pseudo_' num2str(pseudorate(datasetindex)) 'ms']); % ideale classico
    
    %load(['E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\NoNoise\nearFarMP\trackresults_widefilter\trackingResults_data_0_pseudo_' num2str(pseudorate(datasetindex)) 'ms']); %ideale con filtro DLL più largo 20Hz
    %load(['trackingResults_data_ant_1_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
    %load(['trackingResults_data_CN0_100_NoMP_2bit_60s_2_pseudo_' num2str(pseudorate(datasetindex)) 'ms' ]);
    %load(['trackingResults_data_ant_1_thismorning_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
    %load(['trackingResults_data_ant_2_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
    
    %load trackingResults_data_pseudo500ms_16368_41304_noMP_NONOISE_NOSKIP;
    %load 'trackingResults_data_16368_41304_noMP_NONOISE_pseudo_500ms';
    settingsV0=[settingsV0 settings];
    trackResultsV0=[trackResultsV0 trackResults];
    settings.RAIM.enableRAIM=false; settings.enableLAF=false;
    settings.SQIpenalty=1; settings.SQIexclusion=false;
    SQIchannels=zeros(32,10000); settingsMPDD=0;
    LSorWLS=1;
    [navSolutions]=postNavigation(trackResults,settings,settingsMPDD,SQIchannels,LSorWLS);
    %navSolutions= postNavigation(trackResults, settings);
    navSolutionsV0=[navSolutionsV0 navSolutions];
    
    clearvars -except disable navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
    
    if disable==0
        
        load(['trackingResults_data_1_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        %load(['trackingResults_data_ant_2_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        %load(['trackingResults_data_ant_1_noise_AllMP_2bit_60s_5diffWin_pseudo_' num2str(pseudorate(datasetindex)) 'ms' ]);
        %load(['trackingResults_data_ant_3_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        
        
        %load(['trackingResults_data_CN0_100_NoMP_2bit_60s_2_pseudo_' num2str(pseudorate) 'ms' ]);
        %load 'trackingResults_data_16368_41304_CN0_100_NoMP_2bit_60s_pseudo_500ms_elevmask.mat'
        %load(['trackingResults_data_CN0_100_NoMP_2bit_60s_pseudo_' num2str(pseudorate) 'ms' ]);
        %load 'trackingResults_data_16368_41304_noMP_NONOISE_pseudo_500ms_SKIP3s';
        %trackingResults_data_pseudo1000ms_16368_41304_noMP_NONOISE_NOSKIP;
        %trackingResults_data_pseudo1000ms_noise_AllMP_2bit_60s_differentWin
        settingsV1=[settingsV1 settings];
        trackResultsV1=[trackResultsV1 trackResults];
        navSolutions= postNavigation(trackResults, settings);
        navSolutionsV1=[navSolutionsV1 navSolutions];
        
        clearvars -except disable navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
        
        
        load(['trackingResults_data_2_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        
        %load(['trackingResults_data_ant_3_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        %load(['trackingResults_data_NoNoise_MP_2bit_60s_pseudo_' num2str(pseudorate(datasetindex)) 'ms'])
        
        settingsV2=[settingsV2 settings];
        trackResultsV2=[trackResultsV2 trackResults];
        navSolutions= postNavigation(trackResults, settings);
        navSolutionsV2=[navSolutionsV2 navSolutions];
        
        clearvars -except disable navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
        
        %     if disable==0
        
        load(['trackingResults_data_3_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        %     load(['trackingResults_data_ant_1_NOnoise_2MP_2bit_60s_5diffWin_pseudo_' num2str(pseudorate(datasetindex)) 'ms'])
        settingsV3=[settingsV3 settings];
        trackResultsV3=[trackResultsV3 trackResults];
        navSolutions= postNavigation(trackResults, settings);
        navSolutionsV3=[navSolutionsV3 navSolutions];
        
        
        clearvars -except navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
        
        load(['trackingResults_data_4_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        settingsV4=[settingsV4 settings];
        trackResultsV4=[trackResultsV4 trackResults];
        navSolutions= postNavigation(trackResults, settings);
        navSolutionsV4=[navSolutionsV4 navSolutions];
        
        
        
        clearvars -except navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
    end
    
%    load(['/Volumes/Elements/BackUp_Datasets/DataSets/MPMP/newModelMP/5MPtimeWin/differentEph/2bits/NoNoise/trackingResults_maggio2016/trackingResults_data_5_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
     load(['/Volumes/Elements/BackUp_Datasets/DataSets/MPMP/newModelMP/5MPtimeWin/differentEph/2bits/NoNoise/trackresults/trackingResults_data_5_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);

    %load(['trackingResults_data_5_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
    settingsV5=[settingsV5 settings];
    trackResultsV5=[trackResultsV5 trackResults];
    settings.RAIM.enableRAIM=false; settings.enableLAF=false;
    settings.SQIpenalty=1; settings.SQIexclusion=false;
    SQIchannels=zeros(32,10000);
    settingsMPDD=0;
    LSorWLS=1;
    [navSolutions]=postNavigation(trackResults,settings,settingsMPDD,SQIchannels,LSorWLS);
    %navSolutions= postNavigation(trackResults, settings);
    navSolutionsV5=[navSolutionsV5 navSolutions];
    %    end
    
    [firstSubFrame, activeChnList] = findPreambles(trackResultsV0(datasetindex),settingsV0(datasetindex));
    startPreamblems(datasetindex)=firstSubFrame(1);
    
    % test dei singoli datasets
    % plotTracking(1:settings.numberOfChannels, trackResults, settings);
    % [firstSubFrame, activeChnList] = findPreambles(trackResults,settings);
    % startPreamblems=firstSubFrame(1);
    %     navSolutions = postNavigation(trackResults, settings);
    %     [firstSubFrame, activeChnList] = findPreambles(trackResults,settings);
    %     discreteTime=[0:length(navSolutions.channel.rawP(end,:))-1];
    %     xtime=(discreteTime*settings.navSolPeriod+firstSubFrame(1))/1e3;
    %     plotNavigation2(navSolutions, settings,xtime);
    clearvars -except disable navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
    
    
end
%%

%%%%%%%%%%%%%%%%%%%%%
datasetindex=3;
printflag=0;
%% load multipath data
%load ('E:\DataSet\MPMP\newModelMP\5MPtimeWin\params','MPS')
%load('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\NoNoise\nearFarMP\params_2','MPS')
load('/Users/mattiaberardo/Dataset Thuan/2bit/params_2','MPS');


for datasetindex=1:Ndataset
    
    discreteTime=[0:length(navSolutionsV0(datasetindex).channel.rawP(end,:))-1];
    % xtime=(discreteTime+...
    %     floor(settingsV1(datasetindex).msToProcess/settingsV1(datasetindex).navSolPeriod-length(navSolutionsV1(datasetindex).channel.rawP(end,:))))*settingsV1(datasetindex).navSolPeriod/1e3;
    xtime=(discreteTime*settingsV0(datasetindex).navSolPeriod+startPreamblems(datasetindex))/1e3;
    
    %plotNavigation2(navSolutionsV0(datasetindex), settingsV0(datasetindex),xtime,navSolutionsV0(datasetindex));
    %plotNavigation2(navSolutionsV1(datasetindex), settingsV1(datasetindex),xtime,navSolutionsV0(datasetindex));
    %plotNavigation2(navSolutionsV2(datasetindex), settingsV2(datasetindex),xtime,navSolutionsV0(datasetindex));
    %plotNavigation2(navSolutionsV5(datasetindex), settingsV5(datasetindex),xtime,navSolutionsV0(datasetindex));
    
    %cc=trackResultsV5(1:5);
    %dd=settingsV5(1);
    %settingsV5(1).numberOfChannels;
    %plotTracking(1:settingsV5(1).numberOfChannels, cc, dd);
    
    
    
    
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
    %     navSolutionsV5=navSolutionsV2;
    %     settingsV5=settingsV2;
    % xtime=xtime(2:end);
    for index=1:size(navSolutionsV0(datasetindex).channel.PRN,1)
        %[neworder1 I1]=sort(navSolutions1.channel.PRN(:,1));
        %     [neworder2 I2]=sort(navSolutions2.channel.PRN(:,1));
        %     plot(navSolutions2.channel.rawP(I2(index),1:end),'.-')
        
        h=figure;
        subplot(221)
        [neworder1 I1]=sort(navSolutionsV0(datasetindex).channel.PRN(:,1));
        [neworder2 I2]=sort(navSolutionsV5(datasetindex).channel.PRN(:,1));
        %I2=[5 2 4 3 1];%9 15 18 12 14];
        
        plot(xtime,navSolutionsV0(datasetindex).channel.rawP(I1(index),:)-navSolutionsV5(datasetindex).channel.rawP(I2(index),:),'o-',...
            xtime,profileMP1(:)*-36.63,'ro') % *max(abs(navSolutionsV0(datasetindex).channel.rawP(I1(index),:)-navSolutionsV1(datasetindex).channel.rawP(I2(index),1:end)))+1
        grid on
        title(['(pseudorange-pseudorange_{MP}). PRN:' num2str(navSolutionsV5(datasetindex).channel.PRN(I2(index),1))]) % RAW (\rho-\rho_{MP}
        xlabel('Time [s]'), ylabel('\rho-\rho_{MP} [m]')
        
        
        subplot(222)
        plot(xtime,navSolutionsV5(datasetindex).channel.rawP(I2(index),:),'o')
        grid on
        title('Pseudorange') % RAW \rho
        xlabel('Time [s]'), ylabel('\rho [m]')
        
        subplot(224)
        plot(ttime,profileMP2,'ro-')
        grid on
        title('MP profile')
        xlabel('Time [s]'), legend('1 presence of MP/0 absence of MP')
        
        subplot(223)
        plot(xtime,profileMP1(:),'ro-')
        grid on
        title(['MP profile (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms)'])
        xlabel('Time [s]'),legend('1 presence of MP/0 absence of MP')
        
        if printflag==1
            filename=['fig_PRN' num2str(navSolutionsV5(datasetindex).channel.PRN(I2(index),1)) '_pseudorate' num2str(pseudorate(datasetindex))];
            print(h,'-depsc', filename)
        end
    end
    
end




figure, jjj=9; xx=[1:length(trackResultsV5(jjj).I_P)]/1e3;
subplot(411)
plot(xx,trackResultsV5(jjj).I_P,xx,trackResultsV5(jjj).Q_P,'.-'), grid on
title('Correlator V5')
subplot(412)
plot(xx,trackResultsV0(jjj).I_P,xx,trackResultsV0(jjj).Q_P,'.-'), grid on
title('Correlator V0')
subplot(413)
plot(xx,trackResultsV5(jjj).codeFreq-1.023e6,'.-'), grid on
title('Correction V5')
subplot(414)
plot(xx,trackResultsV0(jjj).codeFreq-1.023e6,'.-'), grid on
title('Correction V0')

return
for datasetindex= 1: Ndataset
    
    discreteTime=[0:length(navSolutionsV0(datasetindex).channel.rawP(end,:))-1];
    xtime=discreteTime*settingsV0(datasetindex).navSolPeriod/1e3+startPreamblems(datasetindex)/1e3;
    
    figure
    subplot(321)
    plot(xtime,navSolutionsV0(datasetindex).E-navSolutionsV1(datasetindex).E,'*-',xtime,navSolutionsV0(datasetindex).N-navSolutionsV1(datasetindex).N,'o-',xtime,navSolutionsV0(datasetindex).U-navSolutionsV1(datasetindex).U,'^-')
    title(['Difference in navigation solution (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms). 1 sat with MP'])
    legend('East','North','Up')
    grid on
    
    subplot(322)
    plot(xtime,navSolutionsV0(datasetindex).E-navSolutionsV2(datasetindex).E,'*-',xtime,navSolutionsV0(datasetindex).N-navSolutionsV2(datasetindex).N,'o-',xtime,navSolutionsV0(datasetindex).U-navSolutionsV2(datasetindex).U,'^-')
    title(['Difference in navigation solution (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms). 2 sat with MP'])
    legend('East','North','Up')
    grid on
    
    if disable==0
        
        subplot(323)
        plot(xtime,navSolutionsV0(datasetindex).E-navSolutionsV3(datasetindex).E,'*-',xtime,navSolutionsV0(datasetindex).N-navSolutionsV3(datasetindex).N,'o-',xtime,navSolutionsV0(datasetindex).U-navSolutionsV3(datasetindex).U,'^-')
        title(['Difference in navigation solution (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms). 3 sat with MP'])
        legend('East','North','Up')
        grid on
        
        
        subplot(324)
        plot(xtime,navSolutionsV0(datasetindex).E-navSolutionsV4(datasetindex).E,'*-',xtime,navSolutionsV0(datasetindex).N-navSolutionsV4(datasetindex).N,'o-',xtime,navSolutionsV0(datasetindex).U-navSolutionsV4(datasetindex).U,'^-')
        title(['Difference in navigation solution (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms). 4 sat with MP'])
        legend('East','North','Up')
        grid on
        
        
        %subplot(325)
        plot(xtime,navSolutionsV0(datasetindex).E-navSolutionsV5(datasetindex).E,'*-',xtime,navSolutionsV0(datasetindex).N-navSolutionsV5(datasetindex).N,'o-',xtime,navSolutionsV0(datasetindex).U-navSolutionsV5(datasetindex).U,'^-')%,xtime,profileMP1*80,'r.-','MarkerSize',2)
        title(['Difference in navigation solution (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms). All sat with MP'])
        legend('East','North','Up')
        grid on
    end
    
    
end

timeMP=timeMP/1e3;
fprintf('Multipath istants in seconds:\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n',...
    timeMP(1),timeMP(2),timeMP(2)-timeMP(1),...
    timeMP(3),timeMP(4),timeMP(4)-timeMP(3),...
    timeMP(5),timeMP(6),timeMP(6)-timeMP(5),...
    timeMP(7),timeMP(8),timeMP(8)-timeMP(7),...
    timeMP(9),timeMP(10),timeMP(10)-timeMP(9));
%h = msgbox({'multipath istants in seconds:' num2str(timeMP/1e3) },'MP windows');