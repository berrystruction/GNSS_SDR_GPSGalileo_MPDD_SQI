function [SQIchannels,majorityIndex]=MPDDcomputation(trackResults,settings,settingsMPDD,FLL_results)
SQIchannels=zeros(32,floor(settings.msToProcess/settings.navSolPeriod));

timesize=length(trackResults(1).LAF);

for channelNr = 1:settings.numberOfChannels %[1, 6, 7, 10]
    tic
    %     timesize=length(trackResults(channelNr).LAF);
    seekPRN=settings.seek_sec+FLL_results.seeksec_corr(trackResults(channelNr).PRN);%+ round((20-FLL_results.cnt_skp((trackResults(channelNr).PRN)))*settings.samplingFreq*1e-3);
    timevect=[settingsMPDD.mvgAvgtime:settingsMPDD.mvgAvgtime:timesize*settingsMPDD.mvgAvgtime]/1e3+(seekPRN); % Forse da correggere l asse dei tempi
    wcoeff=zeros(settingsMPDD.M,timesize);
    for timeindex=1:timesize
        wcoeff(:,timeindex)=trackResults(channelNr).LAF(timeindex).w;
    end
    stima_CN0=trackResults(channelNr).CN0.CNo_SNV;
    
    %end
    
    settingsMPDD.MPDDoffset=ceil((((ceil(settings.seek_sec)-settings.seek_sec)/settings.navSolPeriod*1e3)*1e3)/settingsMPDD.mvgAvgtime)+1;
    %settingsMPDD.MPDDoffset=ceil((((ceil(seekPRN)...
    %                        -(seekPRN))/settings.navSolPeriod*1e3)*1e3)/settingsMPDD.mvgAvgtime)+1;
    
    switch settingsMPDD.LAFversion
        case (4)
            alpha=0.4;[indices04, minn04, ~, dictionary, alldistances04]=MPDistanceDetector(alpha,settingsMPDD.M-2,wcoeff(2:end-1,:),timesize,timevect,stima_CN0,settings,settingsMPDD,trackResults(channelNr).PRN);
            alpha=0.5;[indices05, minn05, ~, ~, alldistances05]=MPDistanceDetector(alpha,settingsMPDD.M-2,wcoeff(2:end-1,:),timesize,timevect,stima_CN0,settings,settingsMPDD,trackResults(channelNr).PRN);
            alpha=0.6;[indices06, minn06, ~, ~, alldistances06]=MPDistanceDetector(alpha,settingsMPDD.M-2,wcoeff(2:end-1,:),timesize,timevect,stima_CN0,settings,settingsMPDD,trackResults(channelNr).PRN);
        case (3)
            alpha=0.4;[indices04, minn04, ~, dictionary, alldistances04]=MPDistanceDetector(alpha,settingsMPDD.M,wcoeff(1:end,:),timesize,timevect,stima_CN0,settings,settingsMPDD,trackResults(channelNr).PRN);
            alpha=0.5;[indices05, minn05, ~, ~, alldistances05]=MPDistanceDetector(alpha,settingsMPDD.M,wcoeff(1:end,:),timesize,timevect,stima_CN0,settings,settingsMPDD,trackResults(channelNr).PRN);
            alpha=0.6;[indices06, minn06, ~, ~, alldistances06]=MPDistanceDetector(alpha,settingsMPDD.M,wcoeff(1:end,:),timesize,timevect,stima_CN0,settings,settingsMPDD,trackResults(channelNr).PRN);
    end
    
    % Regole per la detection interDizionario
    if 1==0
        % Ricordarsi di normalizzare così gli indici -> 2*(indices>1)-1
        indices04=2*(indices04>1)-1; indices05=2*(indices05>1)-1; indices06=2*(indices06>1)-1;
        vecIndex =  [indices04; indices05; indices06;];     %vecIndex=[indices03; indices04; indices045; indices05; indices06;];
        majorityIndex=majorityRule(vecIndex,timevect,trackResults(channelNr).PRN);
        DistortionLAFindex=ones(1,timesize);
    else
        %% Altra regola di detection e nuovo parametro dell'indice
        vecIndex =  [indices04; indices05; indices06;];     %vecIndex=[indices03; indices04; indices045; indices05; indices06;];
        minvec   =  [minn04; minn05; minn06];
        LOSdistvec=[alldistances04(1,:); alldistances05(1,:); alldistances06(1,:)]; % vettore con le distanze dalla LOS dei w per ogni dizionario
        
        [absoluteMinDist absoluteMinIndex]=min(minvec); % 1, 2 o 3 nel caso di 3 dizionari
        % Inizializzazione
        absoluteMinCriteria=zeros(1,timesize); DistortionLAFindex=absoluteMinCriteria; LOSdist=DistortionLAFindex;
        
        for timeindex=1:timesize
            absoluteMinCriteria(timeindex)=vecIndex(absoluteMinIndex(timeindex),timeindex);
            % Distanza dalla LOS dopo aver scelto quale vettore dei dizionari usati è +
            % vicino.
            LOSdist(timeindex)=LOSdistvec(absoluteMinIndex(timeindex),timeindex);
            DistortionLAFindex(timeindex)=absoluteMinDist(timeindex)/LOSdist(timeindex); % Se assomiglia al vettore LOS, questo indice è a 1, altrimenti tende a 0
        end
        majorityIndex=2*(absoluteMinCriteria>1)-1;
    end
    
    %% Fine della nuova regola 
    SQIchannels(trackResults(channelNr).PRN,:)=SQIcomputation(trackResults,channelNr,majorityIndex,stima_CN0,settings,settingsMPDD,timesize,seekPRN);
    
    toc
    
    %     % SQI Part
    %     if settingsMPDD.LAFversion==3
    %         [ SQI, judmentfunc ] = MPDD_part2(-majorityIndex,stima_CN0,settings,settingsMPDD,timesize);
    %     else
    %         [ SQI, judmentfunc ] = MPDD_part2(-majorityIndex,stima_CN0,settings,settingsMPDD,timesize,DistortionLAFindex);
    %     end
    %
    %
    %
    %
    %
    %         SQIchannels(trackResults(channelNr).PRN,:)=SQI;
    %         % Print MPDD results
    %         figure
    %
    %         subplot(311)
    %         %     plot(timevect(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1),...
    %         %         majorityIndex(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1),'.','MarkerSize',12),
    %         plot(timevect,majorityIndex,'.','MarkerSize',12),
    %         grid minor
    %         title(['MPDD output for Channel ' num2str(channelNr) ', PRN ' num2str(trackResults(channelNr).PRN)],'interpreter','Latex')
    %         %     title(['Majority Rule for Channel ' num2str(channelNr) ', PRN ' num2str(trackResults(channelNr).PRN)],'interpreter','Latex')
    %         legend('1: MP / -1: NO MP')
    %
    %         subplot(312)
    %         austime=(settingsMPDD.MPDDoffset+[settings.navSolPeriod:settings.navSolPeriod:settings.msToProcess])/1e3+seekPRN; % asse dei tempi
    %         plot(austime,SQI,'ro')
    %         %ylim([-1.5 1.5])
    %         grid minor
    %         title(['SQI for PRN ' num2str(trackResults(channelNr).PRN) ' (with initial offset cutted)'],'interpreter','Latex')
    %         xlabel('Time [s]')
    %         ylim([0 1]);
    %
    %         subplot(313)
    %         plot(timevect,stima_CN0,'.-')
    %         grid on
    %         title('C/N_0')
    %         xlabel('Time [s]')
end







% figure
% %channelNr=6;
% for ttimeindex= 1: length(trackResults(channelNr).LAF) % random time istant index
% subplot(211),
% plot(trackResults(channelNr).LAF(ttimeindex).y,'+-'),  grid minor, hold on, plot(trackResults(channelNr).LAF(ttimeindex).D,'o-'), hold off
% title(['Before MP, Correlation PRN ' num2str(trackResults(channelNr).PRN) ' at time ' num2str(settingsMPDD.mvgAvgtime*1e-3*ttimeindex+settings.seek_sec)])
% subplot(212),
% stem(trackResults(channelNr).LAF(ttimeindex).w),  grid minor
% title('LAF filter Coeff.')
% % subplot(413),
% % plot(Dfinal(:,ttimeindex),'o-'), grid minor, hold on, plot(yFinal(:,ttimeindex),'+-'), hold off
% % subplot(414),
% % stem(wFinal3(:,ttimeindex)),  grid minor
% % title('LAF filter Coeff.')
%
% pause(0.3)
% end

% figure
% subplot(211)
% plot(trackResults(channelNr).carrFreq,'.-')
% subplot(212)
% plot(vec_freq_carr,'.-')


