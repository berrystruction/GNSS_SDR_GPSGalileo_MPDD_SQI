clc
clearvars
close all
addpath('include','geoFunctions');
% addpath('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\2bits\NoNoise\mytrackingresults_widefilter');
%addpath('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\2bits\100dBHz\mytrackingresults_widefilter');
%addpath('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\2bits\60dBHz\mytrackingresults_widefilter');
%addpath('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\2bits\45dBHz\mytrackingresults_widefilter');
addpath('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\7bits\mytrackingresults_widefilter');


% Comparison between data coming from analysis of dataset with different
% pseudorange output rate with/without multipath

% Pseudorate Values available
% 250, 301, 375, 400, 500, 1000, 1044 ms (pseudorange output) 5 sat
%pseudorate=[250 301 375 400 500 1000 1044];
% 103, 301, 500, 1000, 2000 ms (pseudorange output) 4 sat
%pseudorate=[103 301 500 1000 2000];
pseudorate=[1000];
for pseudoindex=1:length(pseudorate)
       
    %Bl=[0.4:0.1:0.6 1  2:2:20 25 30]*10;
    %Bl=[0.5 0.6 1 2:2:10]*10;
    %Bl=[0.4 0.5 0.6 1 2 4 6 8 10 12 18 20 25 30]*10;
    Bl=[2 6 10 20]*10;
    Ndataset=length(Bl);
    
    startPreamblems=zeros(size(Bl));
    
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
    
    
    for datasetindex=1:Ndataset
        
        disp(['loading dataset data #' num2str(datasetindex)]);
        
        %load(['E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\NoNoise\trackresults\trackingResults_data_0_pseudo_' num2str(pseudorate(datasetindex)) 'ms']); % ideale classico
        
        %load(['trackingResults_data_0_pseudo_' num2str(pseudorate(pseudoindex)) 'ms_Bl' num2str(Bl(datasetindex))]); %ideale con filtro DLL più largo 20Hz
        %load(['E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\2bits\100dBHz\mytrackingresults_widefilter\trackingResults_data_0_pseudo_' num2str(pseudorate(pseudoindex)) 'ms_Bl' num2str(Bl(datasetindex))]);
        
        load(['E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\7bits\mytrackingresults_widefilter\trackingResults_data_0_pseudo_' num2str(pseudorate(pseudoindex)) 'ms_Bl' num2str(Bl(datasetindex))]);

        %load(['trackingResults_data_ant_1_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        %load(['trackingResults_data_CN0_100_NoMP_2bit_60s_2_pseudo_' num2str(pseudorate(datasetindex)) 'ms' ]);
        %load(['trackingResults_data_ant_1_thismorning_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        %load(['trackingResults_data_ant_2_pseudo_' num2str(pseudorate(datasetindex)) 'ms']);
        
        %load trackingResults_data_pseudo500ms_16368_41304_noMP_NONOISE_NOSKIP;
        %load 'trackingResults_data_16368_41304_noMP_NONOISE_pseudo_500ms';
        settingsV0=[settingsV0 settings];
        trackResultsV0=[trackResultsV0 trackResults];
        navSolutions= postNavigation(trackResults, settings);
        navSolutionsV0=[navSolutionsV0 navSolutions];
        
        clearvars -except Bl pseudoindex navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
        
        
        load(['trackingResults_data_2_pseudo_' num2str(pseudorate(pseudoindex)) 'ms_Bl' num2str(Bl(datasetindex))]);
        settingsV5=[settingsV5 settings];
        trackResultsV5=[trackResultsV5 trackResults];
        navSolutions= postNavigation(trackResults, settings);
        navSolutionsV5=[navSolutionsV5 navSolutions];
        
        [firstSubFrame, activeChnList] = findPreambles(trackResultsV0(datasetindex),settingsV0(datasetindex));
        startPreamblems(datasetindex)=firstSubFrame(1);
        
        
        clearvars -except Bl pseudoindex navSolutionsV0 settingsV0 trackResultsV0 navSolutionsV1 settingsV1 trackResultsV1 navSolutionsV2 settingsV2 trackResultsV2 navSolutionsV3 settingsV3 trackResultsV3 navSolutionsV4 settingsV4 trackResultsV4 navSolutionsV5 settingsV5 trackResultsV5 pseudorate datasetindex Ndataset startPreamblems;
        
        
    end
    %%
    Bl=Bl/10;
    %%%%%%%%%%%%%%%%%%%%%
    %datasetindex=3;
    printflag=0;
    %% load multipath data
    %load ('E:\DataSet\MPMP\newModelMP\5MPtimeWin\params','MPS')
    load('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\2bits\NoNoise\params_5','MPS')
    %load('E:\DataSet\MPMP\newModelMP\5MPtimeWin\differentEph\NoNoise\nearFarMP\params_2','MPS')
    PRNindex=4;
    trackerror=zeros(length(trackResultsV0(1).I_P),Ndataset);
    for datasetindex=1:Ndataset
        
        discreteTime=[0:length(navSolutionsV0(datasetindex).channel.rawP(end,:))-1];
        xtime=(discreteTime*settingsV0(datasetindex).navSolPeriod+startPreamblems(datasetindex))/1e3;
        %xtime=discreteTime;        
        
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
        for index=1:1%1:size(navSolutionsV0(datasetindex).channel.PRN,1)
            
            h=figure;
            subplot(221)
            [neworder1 I1]=sort(navSolutionsV0(datasetindex).channel.PRN(:,1),'descend');
            [neworder2 I2]=sort(navSolutionsV5(datasetindex).channel.PRN(:,1),'descend');
            %I2=[5 2 4 3 1];%9 15 18 12 14];
            
            diffPseudo=navSolutionsV5(datasetindex).channel.rawP(I2(index),1:end)-navSolutionsV0(datasetindex).channel.rawP(I1(index),:);
            plot(xtime,diffPseudo,'.-',...
                xtime,profileMP1*10,'ro'); % *max(abs(navSolutionsV0(datasetindex).channel.rawP(I1(index),:)-navSolutionsV1(datasetindex).channel.rawP(I2(index),1:end)))+1
            grid on
            title(['(pseudorange-pseudorange_{MP}). PRN ' num2str(navSolutionsV5(datasetindex).channel.PRN(I2(index),1)) '. DLL BW =' num2str(Bl(datasetindex)) 'Hz']) % RAW (\rho-\rho_{MP}
            xlabel('Time [s]')
            ylabel('Diff. between pseudoranges [m]')
            
            
            subplot(222)
            plot(xtime,navSolutionsV5(datasetindex).channel.rawP(I2(index),:),'.')
            grid on
            title('Pseudorange') % RAW \rho
            xlabel('Time [s]')
            ylabel('pseudorange [m]')
            
            subplot(223)
            plot(ttime,profileMP2,'r.-')
            grid on
            title('MP profile')
            xlabel('Time [s]')
            
            subplot(224)
            plot(xtime,profileMP1,'r.-')
            grid on
            title(['MP profile (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms)'])
            xlabel('Time [s]')
            
            if printflag==1
                filename=['fig_PRN' num2str(navSolutionsV5(datasetindex).channel.PRN(I2(index),1)) '_pseudorate' num2str(pseudorate(datasetindex))];
                print(h,'-depsc', filename)
            end
        end
        
        
        
        figure,  xx=[1:length(trackResultsV5(PRNindex).I_P)]/1e3;
        subplot(511)
        plot(xx,trackResultsV5(PRNindex).I_P,xx,trackResultsV5(PRNindex).Q_P,'.-'), grid on
        title(['Prompt Correlators V5, BW = ' num2str(Bl(datasetindex)) 'Hz. PRN ' num2str(trackResultsV5(PRNindex).PRN)])
        xlabel('Time [s]')
        legend('I','Q')
        %plot(xx,trackResultsV0(PRNindex).I_P,xx,trackResultsV0(PRNindex).Q_P,'.-'), grid on
        %title(['Correlators V0, BW = ' num2str(Bl(datasetindex)) 'Hz'])
        subplot(512)
        diffNumSamples=trackResultsV5(PRNindex).absoluteSample-trackResultsV0(PRNindex).absoluteSample;
        %figure
        %plot(xx(1:length(profileMP2)),trackResultsV5(PRNindex).absoluteSample(1:length(profileMP2)),'.-',xx(1:length(profileMP2)),trackResultsV0(PRNindex).absoluteSample(1:length(profileMP2)),'r.-')
        %figure
        plot(xx(1:length(profileMP2)),diffNumSamples(1:length(profileMP2)),'.-',ttime,profileMP2)
        grid on
        title('Difference in # of samples between V5 and V0')
        legend('Difference','MP profile')
        xlabel('Time [s]')
        subplot(513)
        plot(xx,1./trackResultsV0(PRNindex).codeFreq-1./1.023e6,'.-',xx,1./trackResultsV5(PRNindex).codeFreq-1./1.023e6,'.-'), grid on
        %title(['Correction V5, BW = ' num2str(Bl(datasetindex)) 'Hz'])
        %subplot(414)
        %plot(xx,trackResultsV0(PRNindex).codeFreq-1.023e6,'.-'),
        title(['Freq. Code Correction for V0 and V5, BW = ' num2str(Bl(datasetindex)) 'Hz'])
        xlabel('Time [s]')
        ylabel('Freq. Code Corr. [Hz]')
        legend('V0','V5')
        subplot(514)
%         trackV0(:,datasetindex)=trackResultsV0(PRNindex).codeFreq-1.023e6;
%         trackV5(:,datasetindex)=trackResultsV5(PRNindex).codeFreq-1.023e6;
        trackerror(:,datasetindex)=trackResultsV5(PRNindex).codeFreq-1.023e6-(trackResultsV0(PRNindex).codeFreq-1.023e6);
        plot(xx(1:length(profileMP2)),trackerror((1:length(profileMP2)),datasetindex),'.-',ttime,profileMP2), grid on
        title(['Difference between Code Corrections (V5-V0), BW = ' num2str(Bl(datasetindex)) 'Hz'])
        xlabel('Time [s]')
        ylabel('Diff. [Hz]') % Freq. Code Corr.
        legend('Difference','MP profile')
        grid on
        
        subplot(515)
        plot(xtime,diffPseudo,'.-',xtime,profileMP1*10,'ro');
        xlabel('Time [s]')
        title('Diff. [m]')
        title('Diff. between raw Pseudoranges')
        grid on
        PRNindex=PRNindex+5;
        
    end

%% Pseudorange solution comparison   
    figure
    varianz=std(trackerror);
%     varianzV0=var(trackV0);
%     varianzV5=var(trackV5);
%    plot(Bl,varianz,'.-',Bl,varianzV0,'o-',Bl,varianzV5,'^-');
    plot(Bl,varianz,'.-');
    title('Std of the Difference between Code Corrections')
    xlabel('DLL BW [Hz]')
    grid on

%% Navigation solution comparison
    for datasetindex= 1: Ndataset
        
        discreteTime=[0:length(navSolutionsV0(datasetindex).channel.rawP(end,:))-1];
        xtime=discreteTime*settingsV0(datasetindex).navSolPeriod/1e3+startPreamblems(datasetindex)/1e3;
        
        figure
        %subplot(321)
        plot(xtime,navSolutionsV0(datasetindex).X-navSolutionsV5(datasetindex).X,'*-',xtime,navSolutionsV0(datasetindex).Y-navSolutionsV5(datasetindex).Y,'o-',xtime,navSolutionsV0(datasetindex).Z-navSolutionsV5(datasetindex).Z,'^-')
        title(['Difference in Nav. solution (pseudorange rate output: ' num2str(settingsV5(datasetindex).navSolPeriod) ' ms). BW= ' num2str(Bl(datasetindex)) 'Hz.'])
        %legend('East','North','Up')
        legend('X','Y','Z')
        grid on
        
        %varianzE(datasetindex)=var(navSolutionsV0(datasetindex).E-navSolutionsV5(datasetindex).E);
        %varianzN(datasetindex)=var(navSolutionsV0(datasetindex).N-navSolutionsV5(datasetindex).N);
        %varianzU(datasetindex)=var(navSolutionsV0(datasetindex).U-navSolutionsV5(datasetindex).U);
%             varianzE(datasetindex)=mean(((navSolutionsV0(datasetindex).E-navSolutionsV5(datasetindex).E)).^2);
%             varianzN(datasetindex)=mean(((navSolutionsV0(datasetindex).N-navSolutionsV5(datasetindex).N)).^2);
%             varianzU(datasetindex)=mean(((navSolutionsV0(datasetindex).U-navSolutionsV5(datasetindex).U)).^2);
        
    end
    
    
%     figure
%     subplot(311)
%     plot(Bl,varianzE,'.-');
%     title('East variance')
%     xlabel('Hz')
%     grid on
%     subplot(312)
%     plot(Bl,varianzN,'r.-');
%     grid on
%     title('North variance')
%     xlabel('Hz')
%     subplot(313)
%     plot(Bl,varianzU,'.-');
%     title('Up variance')
%     xlabel('Hz')
%     grid on
    
end

timeMP=timeMP/1e3;
fprintf('Multipath istants in seconds:\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n \t [%.3f - %.3f] s, duration: %.3f s\n',...
    timeMP(1),timeMP(2),timeMP(2)-timeMP(1),...
    timeMP(3),timeMP(4),timeMP(4)-timeMP(3),...
    timeMP(5),timeMP(6),timeMP(6)-timeMP(5),...
    timeMP(7),timeMP(8),timeMP(8)-timeMP(7),...
    timeMP(9),timeMP(10),timeMP(10)-timeMP(9));
%h = msgbox({'multipath istants in seconds:' num2str(timeMP/1e3) },'MP windows');
fclose all;