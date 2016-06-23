function [SQIchannels,majorityIndex,LAFerrors,LAFerrorsBiased,Idealerror,ee]=SymmetricGaussianMethod(trackResults,settings,settingsMPDD,FLL_results)
% Gaussian parameters
% Pfa
% Pfa=10^-5;
% GaussThreshold=norminv(Pfa,0,sigma);

timesize=length(trackResults(1).LAF);
SQIchannels=zeros(32,floor(settings.msToProcess/settings.navSolPeriod));

LAFerrors=zeros(settings.numberOfChannels,timesize);
LAFerrorsBiased=LAFerrors;
Idealerror=LAFerrors;
len=length(trackResults(1).LAF(1).e(1:end));
NpointsLR=(len-1)/2;
len=length(trackResults(1).LAF(1).e(ceil(len/2)-NpointsLR:ceil(len/2)+NpointsLR));
dd=zeros((len-1)*timesize,1);
ee=zeros(settings.numberOfChannels,length(dd));
%cc=zeros(len,1);

for channelNr = 1:settings.numberOfChannels
    tic
    majorityIndex=zeros(1,timesize);
    %stima_CN0=trackResults(channelNr).CN0.CNo_SNV;
    %     timesize=length(trackResults(channelNr).LAF);
    indice=1;
    %cc=zeros(len,1);
    for timeindex=1:timesize
        %approximErrLAF=trackResults(channelNr).LAF(timeindex).e(9:33);
        %cc=[approximErrLAF([1:floor(len/2)]); -approximErrLAF([floor(len/2)+2:end])];
        
        approximErrLAF=trackResults(channelNr).LAF(timeindex).e(ceil(len/2)-NpointsLR:ceil(len/2)+NpointsLR);
        cc=[approximErrLAF(ceil(len/2)-NpointsLR:ceil(len/2)-1); -approximErrLAF(ceil(len/2)+1:ceil(len/2)+NpointsLR)];
        

        LAFerrors(channelNr,timeindex)=mean(cc(:));
        
        cc_biased=[approximErrLAF([1:floor(len/2)]); approximErrLAF(ceil(len/2)); -approximErrLAF([floor(len/2)+2:end])];
        LAFerrorsBiased(channelNr,timeindex)=mean(cc_biased);
        
        approximErr=trackResults(channelNr).LAF(timeindex).D-trackResults(channelNr).LAF(timeindex).U(:,ceil(settingsMPDD.M/2));
        Idealerror(channelNr,timeindex)=mean([approximErr([1:floor(len/2)]); -approximErr([floor(len/2)+2:end])]);
        
        dd(indice:indice+len-2)= cc;
        indice=indice+len-1;
        %sigma=std(LAFerrors(timeindex-settingsMPDD.mvgAvgtime+1:timeindex));
        %         sigma2=.001;%var(LAFerrors(1,1:20));
        %         if timeindex>10
        %             %T_x=mean(LAFerrors(channelNr,timeindex-10+1:timeindex));
        %         sigma2=var([cc(:,1); cc(:,2); cc(:,3);cc(:,4);cc(:,5);cc(:,6);cc(:,7);cc(:,8);cc(:,9);cc(:,10)]);
        %         T_x=LAFerrors(channelNr,timeindex);
        %
        %         else
        %             T_x=0;
        %         end
        % %         if stima_CN0(timeindex)>40 % dBHz
        %               Pfa=10^-15;
        % %             sigma=stima_CN0(timeindex);
        % %         else
        % %             Pfa=10^-5;
        % %             sigma=stima_CN0(timeindex);
        % %         end
        %         GaussThreshold=sqrt(sigma2/(len-1))*qfuncinv(Pfa);
        %         majorityIndex(timeindex)=2*(T_x/sqrt(sigma2)>GaussThreshold)-1;
    end
    
    ee(channelNr,:)=dd';
    
    seekPRN=settings.seek_sec+FLL_results.seeksec_corr(trackResults(channelNr).PRN);%+ round((20-FLL_results.cnt_skp((trackResults(channelNr).PRN)))*settings.samplingFreq*1e-3);
    settingsMPDD.MPDDoffset=ceil((((ceil(settings.seek_sec)-settings.seek_sec)/settings.navSolPeriod*1e3)*1e3)/settingsMPDD.mvgAvgtime)+1;
    %settingsMPDD.MPDDoffset=ceil((((ceil(seekPRN)...
    %                        -(seekPRN))/settings.navSolPeriod*1e3)*1e3)/settingsMPDD.mvgAvgtime)+1;
    
    %% Fine della nuova regola
    %SQIchannels=SQIcomputation(trackResults,channelNr,majorityIndex,stima_CN0,settings,settingsMPDD,timesize,seekPRN);
    
    toc
    
end


