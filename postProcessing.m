% Script postProcessing.m processes the raw signal from the specified data
% file (in settings) operating on blocks of 37 seconds of data.
%
% First it runs acquisition code identifying the satellites in the file,
% then the code and carrier for each of the satellites are tracked, storing
% the 1msec accumulations.  After processing all satellites in the 37 sec
% data block, then postNavigation is called. It calculates pseudoranges
% and attempts a position solutions. At the end plots are made for that
% block of data.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis, Dennis M. Akos
% Some ideas by Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%                         THE SCRIPT "RECIPE"
%
% The purpose of this script is to combine all parts of the software
% receiver.
%
% 1.1) Open the data file for the processing and seek to desired point.
%
% 2.1) Acquire satellites
%
% 3.1) Initialize channels (preRun.m).
% 3.2) Pass the channel structure and the file identifier to the tracking
% function. It will read and process the data. The tracking results are
% stored in the trackResults structure. The results can be accessed this
% way (the results are stored each millisecond):
% trackResults(channelNumber).XXX(fromMillisecond : toMillisecond), where
% XXX is a field name of the result (e.g. I_P, codePhase etc.)
%
% 4) Pass tracking results to the navigation solution function. It will
% decode navigation messages, find satellite positions, measure
% pseudoranges and find receiver position.
%
% 5) Plot the results.

%% Initialization =========================================================
disp ('Starting processing...');

[fid, message] = fopen(fileName, 'rb');

%If success, then process the data
if (fid > 0)
    
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    status=fseek(fid, settings.skipNumberOfBytes, 'bof');
    if (status)
        error('fseek error')
    end
    
    %% Acquisition ============================================================
    
    % Do acquisition if it is not disabled in settings or if the variable
    % acqResults does not exist.
    if ((settings.skipAcquisition == 0) || ~exist('acqResults', 'var'))
        
        % Find number of samples per spreading code
        samplesPerCode = round(settings.samplingFreq / ...
            (settings.codeFreqBasis / settings.codeLength));
        
        % Read data for acquisition. 11ms of signal are needed for the fine
        % frequency estimation
        %         data = fread(fid, 11*samplesPerCode, settings.dataType)';
        
        acqNsamp = 1*samplesPerCode*(settings.Non_Coh_Sums)*settings.NSample;
        [data] = fread(fid,acqNsamp, settings.dataType)';
        
        %if (strcmp(settings.signalType,'IQ')==1) % IQ sampling
        %    [data] = fread(fid,acqNsamp,settings.dataType)'; % read 4 times the number of samples for the baseband complex (16bitx2) signal
        %else % IF sampling
            %acqNsamp = 1*samplesPerCode*(settings.Non_Coh_Sums)*settings.NSample;
        %    [data] = fread(fid,acqNsamp, settings.dataType)';
        %end
        
        if (length(data)~=acqNsamp)
            disp('Not able to read the specified number of samples, exiting...')
            fclose(file);
            return
        end
        
        
        
        %--- Do the acquisition -------------------------------------------
        disp ('   Acquiring satellites...');
        %acqResults = acquisition(data, settings);
        acqResults=myacquisition(data,settings);
        
        plotAcquisition(acqResults,settings);
        
    end
    
    
    %     %% FLL ====================================================================
    %    FLL_results=trackcarrFLL(fid,acqResults,settings);
    
    % DA SISTEMARE NON ANCORA CORRETTO PER BENE
    % settings.seek_sec=settings.seek_sec+(FLL_results(settings.acqSatelliteList(1)).phaseFLL/settings.samplingFreq-settings.seek_sec); % correction of skipped seconds
    
    
    %% Initialize channels and prepare for the run ============================
    
    % Start further processing only if a GNSS signal was acquired (the
    % field FREQUENCY will be set to 0 for all not acquired signals)
    %    if (any(acqResults.carrFreq))
    if (any(acqResults.peakMetric >= settings.acqThreshold)) %%% Ottimo sia per BS che per Realsamples
        [channel,FLL_results] = preRun(fid,acqResults,settings); % ... And integration of FLL
        showChannelStatus(channel,settings);
    else
        % No satellites to track, exit
        disp('No GNSS signals detected, signal processing finished.');
        trackResults = [];
        return;
    end
    
    
    % ====================================================================
    %% Track the signal =======================================================
    startTime = now;
    disp (['   Tracking started at ', datestr(startTime)]);
    
    %settings.codeLength         = 1023*5;
    
    % Process all channels for given data block
    [trackResults, channel] = tracking(fid,channel,settings,settingsMPDD,FLL_results);
    
    % Close the data file
    fclose(fid);
    
    disp(['   Tracking is over (elapsed time ', ...
        datestr(now - startTime, 13), ')'])
    
    
    %% MPDD computation in order to assit RAIM ================================
    if settings.enableLAF==true %&& 1==0
        
        %%trackResults2=trackResults;
        %%trackResults=trackResults2([1:8 10]); 
        %%settings.numberOfChannels=length(trackResults);
        %[SQIchannels, MPDDpoints]=MPDDcomputation(trackResults,settings,settingsMPDD,FLL_results);
        
%         [SQIchannels,MPDDpoints,LAFerrors,LAFerrorsBiased,Idealerror]=SymmetricGaussianMethod(trackResults,settings,settingsMPDD,FLL_results);
%         for chn=1:length(trackResults), figure, plot(LAFerrors(chn,:),'bo-'), hold on, plot(Idealerror(chn,:),'r.-'), hold off,grid minor,title(['Channel ' num2str(chn)]), legend('Unbiased error','Biased error'), end
%         [h,p] = chi2gof(LAFerrors,'Alpha',0.05); h
%         [h,p] = chi2gof(LAFerrorsBiased,'Alpha',0.05); h
       %% 
       close all
       significanceLevel=10/100;
       [SQIchannels,Idealerror,MPresults,delayestimation]=iterativeSGM(trackResults,settings,settingsMPDD,FLL_results,significanceLevel);
       %%
%        figure, chn=1; lenw=length(trackResults(chn).LAF(1).w); lendiff=length(trackResults(chn).LAF(1).D)-length(trackResults(chn).LAF(1).d);
%        for kk=1:5:length(trackResults(1).LAF)%/2+20, 
%            subplot(221), plot(trackResults(chn).LAF(kk).D,'o-'), hold on, plot(trackResults(chn).LAF(kk).y,'*-'), hold off
%            title(['Time ' num2str(kk*settingsMPDD.mvgAvgtime*1e-3+settings.seek_sec) ' PRN ' num2str(trackResults(chn).PRN)]), grid on,
%            subplot(222), 
%            ccU=trackResults(chn).LAF(kk).U(:,[1:floor(lenw/2) floor(lenw/2)+2:end]);
%            ccw=trackResults(chn).LAF(kk).w([1:floor(lenw/2) floor(lenw/2)+2:end]);
%            ccUU=trackResults(chn).LAF(kk).U(:,floor(lenw/2)+1)*trackResults(chn).LAF(kk).w(floor(lenw/2)+1);
%            plot(ccU*ccw,'o-'), grid on, hold on, plot(ccUU,'*-'), grid on, hold on, plot(ccU*ccw+ccUU,'+-'), grid on, legend('Residual','U_{los}*w','y'), hold off
%            %for kk2=1:lenw, plot(trackResults(chn).LAF(kk).U(:,kk2)*trackResults(chn).LAF(kk).w(kk2),'o-'), hold on, grid on, end, hold off
%            subplot(223), stem(trackResults(chn).LAF(kk).w), grid on,  
%            subplot(224), plot(trackResults(chn).LAF(kk).D(lendiff/2+1:end-lendiff/2+1)-trackResults(chn).LAF(kk).y(lendiff/2+1:end-lendiff/2+1),'o-'), grid on
%            pause(.3)
%        end
       
    else
        SQIchannels=zeros(32,floor(settings.msToProcess/settings.navSolPeriod));
    end
    
    %% Auto save the acquisition & tracking results to a file to allow
    % running the positioning solution afterwards.
    
    [tok, remain]=strtok(fileName(end-4:-1:1),'\/'); tok=tok(end:-1:1);
    myfolder=remain(end:-1:1);
    %myfolder=[remain(end:-1:1) 'mytrackingResults'];
    %mkdir(myfolder);
    
    %myfolder=[myfolder num2str(round(rem(now,1)*1e4)) '_trackingResults_' tok '_pseudo_' num2str(settings.navSolPeriod) 'ms'];
    % New file name template
    myfolder=[myfolder num2str(round(rem(now,1)*1e4)) '_trackingRes_' tok '_' num2str(1/floor(settings.navSolPeriod/1000)) 'Hz_' num2str(settings.msToProcess)  'ms'];
    %_Bl' num2str(settings.dllNoiseBandwidth*10)];
    save(myfolder)
    %, ...'trackResults', 'settings','settingsMPDD', 'acqResults', 'channel','FLL_results','SQIchannels');
    
    disp('   Saving Acq & Tracking results to .mat file ') % "trackingResults.mat"
    
    %% Calculate navigation solutions =========================================
    disp('   Calculating navigation solutions...');
    %%
    initConstant
    settings.useIonoCorr=0;
    %settings.DOPcontrol=false;
    %settings.enableLAF=false;
    
    settings.enableSQI=false; % To handle unreliable solutions
    settings.SQIpenalty=0.0; % SQI covariance penalty factor
    settings.SQIexclusion=true; % SQI controlled exclusion
    
    %settings.RAIM.enableRAIM=true;
    %settings.RAIM.type=0; % 0-RS (Residual method), 1-SS (Solution separation method)
    %settings.RAIM.RecomputePVT=5;
    
    kalmanSettings.staticposition=false;
    
    trackResultsNsat=trackResults;%([1 6 3 2 9 10]);%([1:8 9 10]);%([10 6 9 3 2 1]);%([1 3 4 10 7]);%([1 2 6 8 9]);%([1:7]);%([1:6  8:10]);%([1:3 5 8 9]);%([1:6 8]);%([1:8]);%([ 1:7 9 10]);%([1 7 8 9 10 11 12 13]);
    
%     for chn=[1:length(trackResultsNsat)]
%         figure, plot(abs(trackResultsNsat(chn).I_P),'b'), hold on,
%         plot(abs(trackResultsNsat(chn).I_E),'r'), hold on,...
%             plot(abs(trackResultsNsat(chn).I_L),'g'), grid on, title(['channel ' num2str(chn) ' PRN ' num2str(trackResultsNsat(chn).PRN)]), legend('Prompt','Early','Late');
%         %plot(trackResultsNsat(chn).CN0.CNo_SNV,'.-') , hold on
%     end
%     %% Code freq. corrections plot
%     close all, clc, chn=2;
%     xtime=[1:length(trackResults(chn).codeFreq)]*1/settings.codeFreqBasis*settings.codeLength;%+settings.seek_sec;
%     figure, subplot(211), plot(xtime,-(1/settings.codeFreqBasis - 1./trackResults(chn).codeFreq),xtime,mean(-(1/settings.codeFreqBasis - 1./trackResults(chn).codeFreq))*ones(size(xtime)),'r-'), xlabel('Time [s]'), ylabel('Corrections [s]')
%     grid on, title('T_{chip} corrections'), subplot(212), plot(xtime,-(trackResults(chn).codeFreq - settings.codeFreqBasis),xtime,mean(-(trackResults(chn).codeFreq - settings.codeFreqBasis))*ones(size(xtime)),'r'), grid on, xlabel('Time [s]'), ylabel('Corrections [Hz]'), title('NCO corrections')
%     
    LSorWLS=0;  % parameters: 1 -> LS, 0 -> WLS
    
    for pvt2=1:1
        
        if pvt2==1 % W/LS
            [navSolutions,~,RAIMresults,unreliableSol,TGlobalTest,TlocalTest,GlobalThres,remainingSV,SVexcluded]=...
                postNavigation(trackResultsNsat,settings,settingsMPDD,SQIchannels,LSorWLS);
            
            %kalmanSettings.startKalman=Inf;
            %[navSolutions,RAIMresults,unreliableSol,TGlobalTest,TlocalTest,GlobalThres,remainingSV,SVexcluded]=PVT_Kalman(trackResultsNsat,settings,kalmanSettings,LSorWLS,SQIchannels,settingsMPDD.SQIthreshold);
            
            plotNavigation(navSolutions,settings,myfolder);
            plotNavigationNoref(navSolutions);
            
           %%
            if isnan(settings.truePosition.E)==0 % If I have a reference position
                figure
                ENUerrors(1,:)=(navSolutions.E - settings.truePosition.E);
                ENUerrors(2,:)=(navSolutions.N - settings.truePosition.N);
                ENUerrors(3,:)=(navSolutions.U - settings.truePosition.U);
                [aa,bb]=hist(ENUerrors(3,:),30);
                bar(bb,aa), grid on, title('Histogram of the vertical error'), xlabel('error [m]')
                figure
                [cc,dd]=hist(sqrt(ENUerrors(1,:).^2+ENUerrors(2,:).^2),30); bar(dd,cc), grid on, title('Histogram of the horizontal error'), xlabel('error [m]')
            end
%             figure
%             subplot(211)
%             plot(sqrt(ENUerrors(1,:).^2+ENUerrors(2,:).^2),[navSolutions.PL(:).HPL]), grid on
%             subplot(212)
%             plot(sqrt(ENUerrors(1,:).^2+ENUerrors(2,:).^2),[navSolutions.PL(:).HPL2]), grid on
            %%
            %figure,
            %ausiliarVar=[RAIMresults([1:end]).w];
            %for kk=1:length(trackResultsNsat), plot(ausiliarVar(kk,:),'.-'), grid on, hold on, end, hold off ,legend('1','2','3','4','5','6','7','8','9') ,title('LS residuals')
            % Plot comparisons
            % navSolComparison(navSolutions,upos,settings.navSolPeriod);
        else % KF
            kalmanSettings.startKalman=7;
            [navSolutionsKF,RAIMresultsKF,unreliableSol,TGlobalTest,TlocalTest,GlobalThres,remainingSV,SVexcluded]=...
                PVT_Kalman(trackResultsNsat,settings,kalmanSettings,LSorWLS,SQIchannels,settingsMPDD.SQIthreshold);
            %PVT_Kalman_FallBack(trackResultsNsat,settings,kalmanSettings,LSorWLS,SQIchannels,settingsMPDD.SQIthreshold);
            % fallback for the work with Hieu
            plotNavigation(navSolutionsKF,settings,myfolder);
            plotNavigationNoref(navSolutionsKF);
            
            %navSolutions=navSolutionsKF;
            
            %figure,
            %ausiliarVar=[RAIMresultsKF([1:end]).w];
            %for kk=1:length(trackResultsNsat), plot(ausiliarVar(kk,:),'.-'), grid on, hold on, end, hold off, legend('1','2','3','4','5','6','7','8','9') ,title('Kalman residuals')
            % Plot comparisons
            %navSolComparison(navSolutionsKF,upos,settings.navSolPeriod,navSolutions);
        end
        
        disp('PRNs:')
        [trackResultsNsat(:).PRN]
        
        % RAIM
        if settings.numberOfChannels>=6 && settings.RAIM.enableRAIM==true && exist('navSolutions','var') %(isfield('navSolutions','X') || isfield('navSolutionsKF','X') )% and NavSolutions variable exists
            RAIM_Execution(TGlobalTest,TlocalTest,GlobalThres,unreliableSol,settings,navSolutions)
        end
        disp('   Processing is complete for this data block');
        
        % Plot all results ===================================================
        disp ('   Ploting results...');
        if settings.plotTracking
            plotTracking(1:settings.numberOfChannels, trackResults, settings);
        end
        
        % Plotting # of satellite not excluded vs time
        if isempty(navSolutions)==0
            figure
            xtime=([0:length(navSolutions.E(end,:))-1]+floor(settings.msToProcess/settings.navSolPeriod-length(navSolutions.E(end,:))))*settings.navSolPeriod/1e3;
            plot(xtime,remainingSV,'o'),grid minor,xlabel('Time [s]'),title('Remaining SV per epoch')
            
            % Residual pseudorange trend
            if settings.RAIM.enableRAIM==true && isempty(RAIMresults)==0 && 1==0
                %%
                figure
                for PRNindex=[1:1]%length(trackResultsNsat)]
                    residualsinglePRN=zeros(length(RAIMresults(1).w(:)),length(RAIMresults));
                    %                         istantaneouscov=zeros(length(RAIMresults(1).w(:)),length(RAIMresults(1).w(:)),length(RAIMresults));
                    indexx2=1;
                    timewin=5;
                    for indexx=1:length(RAIMresults)
                        %if exist('RAIMresults(indexx).w(PRNindex)','var')==1
                        residualsinglePRN(:,indexx)=RAIMresults(indexx).w(:);
                        %                             istantaneouscov(:,:,indexx2)=istantaneouscov(:,:,indexx2)+residualsinglePRN(:,indexx)*residualsinglePRN(:,indexx)';
                        %                             if rem(indexx,timewin)==0
                        %                                 istantaneouscov(:,:,indexx2)=istantaneouscov(:,:,indexx2)/timewin;
                        %                                 indexx2=indexx2+1;
                        %                                 [W,L,expl] = pcacov(istantaneouscov(:,:,1));
                        %                                 %[W,score,L]=pca(residualsinglePRN'); L/sum(L)*100
                        %                                 %biplot(W(:,1:2),'Scores',score(:,1:2))
                        %                             end
                        
                        %else
                        %   residualsinglePRN(indexx)=NaN;
                        %end
                    end
                    %simple = tsmovavg((residualsinglePRN),'s',10,1);
                    %figure
                    %detrend_sdata = detrend(residualsinglePRN);
                    
                    %Z=residualsinglePRN'*W(:,1:2);
                    %figure
                    plot(xtime,residualsinglePRN(:,:),'.-')%,xtime,simple,'o-'),
                    %ylim([-50 25]), hold on,
                end
                legend(strsplit([num2str([trackResultsNsat(:).PRN])]))
                grid minor, title('PRN Residual pseudorange trend'), xlabel('Time [s]'), ylabel('Amplitude [m]')
                %
                figure,
                rawPdiff=navSolutions.channel.rawP(:,3:end)-navSolutions.channel.rawP(:,2:end-1);
                detrend_sdata = detrend(rawPdiff')'; % Remove linear trends
                subplot(211)
                plot(xtime(3:end),rawPdiff,'.-'), grid on, title('Residual pseudorange in time'), legend(strsplit([num2str([navSolutions.channel.PRN(:,1)'])])),xlabel('Time [s]'),ylabel('Amplitude [m]')
                subplot(212), plot(xtime(3:end),detrend_sdata,'.-'), grid on, title('Detrended Residual pseudorange in time'),xlabel('Time [s]'),
                
                figure,
                for kk=1:length([trackResultsNsat(:).PRN])
                    PRN=trackResultsNsat(kk).PRN;
                    plot(xtime,SQIchannels(PRN,6:end-2),'o-'), grid minor, hold on
                end
                legend(strsplit([num2str([trackResultsNsat(:).PRN])]))
                
            end
            %% Variance estimation for protection level computation
            % % % % %             if settings.enableLAF==true
            % % % % %                 varianceMtx = PseudoVarianceComputation(settings,settingsMPDD.SQIthreshold,navSolutions.channel.correctedP,SQIchannels);
            % % % % %                 figure,
            % % % % %                 plot(xtime,varianceMtx(1,:))
            % % % % %                 grid on
            % % % % %             end
            
            %figure,
            
            %         for kk=9:9%[1:2 4 5 7 8] % Alternative MP detection method??
            %             init=1;
            %             finit=100;
            %             for init=50:100:900
            %                 init
            %                 nn=(sqrt(sum([trackResults(kk).LAF(init:finit).e].^2,1)));
            %                 nn=nn-mean(nn); mean(nn)
            %                 stdnn=std(nn);
            %                 figure
            %                 subplot(211)
            %                 plot(nn), grid on,legend(['PRN ' num2str(trackResults(kk).PRN)]), title('norm of approximation Error vector for the LAF')
            %                 subplot(212)
            %                 nbins=30;
            %                 [counts,xcent]=hist(nn,nbins);
            %                 bar(xcent,counts/trapz(xcent,counts)), grid minor, %hold on,
            %                 [~,~,flagg]=gaussianity(nn,nbins);
            %                 %g=1/sqrt(2*pi)*exp(-0.5*xcent.^2);%# pdf of the normal distribution
            %                 %plot(xcent,g,'r-'); hold off;
            %                 %init=finit+1;
            %                 finit=finit+100;
            %             end
            %
            %         end
            
        end
        
        if isempty(SVexcluded)==0
            SVexcluded(1:end).PRN
        else
            disp('No satellites excluded')
        end
        
        disp('-------------------------------------------------------------')
        
    end
    
    disp('Post processing of the signal is over.');
    
else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', fileName, message);
end % if (fid > 0)
