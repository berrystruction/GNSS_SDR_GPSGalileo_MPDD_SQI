function settingsMPDD = InitsettingsMPDD(settings)
%% Processing settingsMPDD ====================================================
% version can be:
% 2- normal LAF 
% 3 CLS- autocorr. LAF
% 4 LS- autocorr. LAF 
% 5 IHT 
% 6 LS con equality constraint
version=6; 
settingsMPDD.LAFversion=version;
%% MULTICORRELATOR AND LAF ====================================================
settingsMPDD.mvgAvgtime=500; % 150; % finestra di moving average; if mvgAvgtime==1 attivo anche il postavgLAF
settingsMPDD.multicorr.maxlag=1.0;%0.7;%1.7;%0.99; % in terms of Tchip
%Npoints=15; %16 � la risoluzione del Ts, (se aumento devo allungare il filtro) %lags for the LAF correlations
% calcolo il Npoints minimo x avere risoluzione maggiore di Ts, posso
% calcolarla in 2 modi
settingsMPDD.multicorr.Tres_Tsamp=0.937;%0.979;%1.373;%3.5;%2 %1.47;%0.98;%0.97 % percent. Resolution time vs sampling time
%settingsMPDD.multicorr.Npoints=round(settingsMPDD.multicorr.maxlag/1.023e6/(settingsMPDD.multicorr.Tres_Tsamp/settings.samplingFreq));
settingsMPDD.multicorr.Npoints=9;%9;
Npoints=settingsMPDD.multicorr.Npoints;

%Npoints=ceil(maxlag/1.023e6*f_sampling); % modo 1 risoluzione superiore a quella del campionamento
%(2*maxlag*1/1.023e6/(1/f_sampling*2)-0.5) % modo 2 risoluzione superiore a quella del campionamento
%Npoints=floor(maxlag/1.023e6*f_sampling)-2; % risoluzione inferiore a quella del campionamento


% M=round((2*Npoints+1)*23/100); % # of tap for the LAF filter
% M_Nratio=M/(2*Npoints+1);

%% multicorrelator resolution, enable this part in order to monitor the resolution changing
settingsMPDD.multicorr.multicorr_resolution=settingsMPDD.multicorr.maxlag*1/settings.codeFreqBasis/(Npoints); % multicorr_resolution=2*maxlag*1/1.023e6/(2*Npoints-1)*1e6
% M_visibitlityRange=multicorr_resolution*M*1e6;
% oppure imposto il visibility range e scelgo M di conseguenza
settingsMPDD.multicorr.M_visibitlityRange=settingsMPDD.multicorr.maxlag;%3.0; % voluto
settingsMPDD.M=ceil(settingsMPDD.multicorr.M_visibitlityRange/1.023e6/settingsMPDD.multicorr.multicorr_resolution);

%settingsMPDD.M=21;% 17

settingsMPDD.M_Nratio=settingsMPDD.M/(2*Npoints+1);

%soglia=.27;
%if settingsMPDD.M_Nratio>=soglia
%    fprintf('\n M=%d, maybe is too big respect to N! MNR=%f!\n',settingsMPDD.M,settingsMPDD.M_Nratio);
%     while settingsMPDD.M_Nratio>=soglia
%         settingsMPDD.M=settingsMPDD.M-1;
%         settingsMPDD.M_Nratio=settingsMPDD.M/(2*Npoints+1);
%     end
%end
settingsMPDD.multicorr.M_visibitlityRange=settingsMPDD.multicorr.multicorr_resolution*(settingsMPDD.M-2)*1e6; % ricalcolato dopo aver impostato M
settingsMPDD.multicorr.upperBoundVR=(2*Npoints+1)*0.5*settingsMPDD.multicorr.multicorr_resolution*settings.codeFreqBasis; % max VR, VR < upperbound
%%

%%
settingsMPDD.mvgAvgtimeflag=false; % flag Avgpost-LAF
if settingsMPDD.mvgAvgtimeflag==true
    settingsMPDD.mvgAvgtime2=2; % window post LAF
else
    settingsMPDD.mvgAvgtime2=1;
end
settingsMPDD.datasaveFlag=true; % to enable save

%% print
if settings.enableLAF==true
    fprintf('Correlation spacing is: %1.2f with:\n',2*settingsMPDD.multicorr.maxlag)
    fprintf(' 2Npoints+1=%d,\n 2Npoints+1-Npoints/2=%d,\n 2Npoints+1-Npoints/2-M+1=%d\n',2*Npoints+1,2*Npoints+1-floor(Npoints/2),2*Npoints+1-floor(Npoints/2)-settingsMPDD.M+1)
    fprintf('\n Npoints=%d, \t\t M=%d, \t\t  MNR=%f, \n resolution of multicorrelator=%f us,\n multicorr_resolution/sampling ratio=%f,\n',...
        Npoints,settingsMPDD.M,settingsMPDD.M_Nratio,settingsMPDD.multicorr.multicorr_resolution*1e6,settingsMPDD.multicorr.multicorr_resolution*settings.samplingFreq);
    fprintf(' Visibility_range (M-2)=%1.2f Chip,\t Upperbound VR=%f,\n',settingsMPDD.multicorr.M_visibitlityRange,settingsMPDD.multicorr.upperBoundVR);
    fprintf(' Visibility_range (M)=%1.2f Chip,\n',settingsMPDD.multicorr.M_visibitlityRange/(settingsMPDD.M-2)*settingsMPDD.M);
    fprintf(' Upperbound M=%3.2f\n',settingsMPDD.multicorr.upperBoundVR/settings.codeFreqBasis/settingsMPDD.multicorr.multicorr_resolution);
    fprintf(' LAF version=%d\n\n',settingsMPDD.LAFversion);
end

%% Acquisition settingsMPDD ===================================================
% % Skips acquisition in the script postProcessing.m if set to 1
% settingsMPDD.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
%settingsMPDD.acqSatelliteList   = [1:32];%1:32;         %[PRN numbers]
% Band around IF to search for satellite signal. Depends on max Doppler
%settingsMPDD.acqSearchBand      = 14;           %[kHz]
% Threshold for the signal presence decision rule
%settingsMPDD.acqThreshold       = 2.0;

%% Tracking loops settingsMPDD ================================================
% Code tracking loop parameters
% settingsMPDD.dllDampingRatio         = 0.7; %1/sqrt(2);
% settingsMPDD.dllNoiseBandwidth       = 2;       %[Hz]
% settingsMPDD.dllCorrelatorSpacing    = 0.5;     %[chips]
%
% % Carrier tracking loop parameters
% settingsMPDD.pllDampingRatio         = 0.7; %1/sqrt(2);
% settingsMPDD.pllNoiseBandwidth       = 25;      %[Hz]
%% Navigation solution settingsMPDD ===========================================
% Period for calculating pseudoranges and position
%settingsMPDD.navSolPeriod       = 1000;         %[ms]
settingsMPDD.MPDDNpoints=settings.navSolPeriod/settingsMPDD.mvgAvgtime;
%% CN0 estimation parameters ==============================================
settingsMPDD.CNo_WINDOW=settingsMPDD.mvgAvgtime;%settingsMPDD.mvgAvgtime;

%% parameter MPDD part 2
settingsMPDD.minCN0considered=40; %[dBHz]
settingsMPDD.MPDDoffset=NaN; % [MPDD point] about 500 ms (offset to put istant of PVT publication in the middle of the MPDD epoch)
settingsMPDD.plotMPDD=0;
settingsMPDD.SQIthreshold=0.65;
% % Elevation mask to exclude signals from satellites at low elevation
% settingsMPDD.elevationMask      = 10;           %[degrees 0 - 90]
% % Enable/dissable use of tropospheric correction
% settingsMPDD.useTropCorr        = 1;            % 0 - Off
%                                             % 1 - On
%
% % True position of the antenna in UTM system (if known). Otherwise enter
% % all NaN's and mean position will be used as a reference .
% settingsMPDD.truePosition.E     = nan;
% settingsMPDD.truePosition.N     = nan;
% settingsMPDD.truePosition.U     = nan;
%
% %% Plot settingsMPDD ==========================================================
% % Enable/disable plotting of the tracking results for each channel
% settingsMPDD.plotTracking       = 1;            % 0 - Off
%                                             % 1 - On
%
% %% Constants ==============================================================
%
%settingsMPDD.c                  = 299792458;    % The speed of light, [m/s]
%settingsMPDD.ProcessingDate     = NaN;
% settingsMPDD.startOffset        = 68.802;       %[ms] Initial sign. travel time
%settingsMPDD.Nbits              = NaN;
end

