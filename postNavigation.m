function [navSolutions, eph, RAIMresults, unreliableSol, TGlobalTest, TlocalTest, GlobalThres, remainingSV, SVexcluded, IONOdelay] = postNavigation(trackResults,settings,settingsMPDD,SQIchannels,LSorWLS)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis with help from Kristin Larson
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

%CVS record:
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $

%% Check is there enough data to obtain any navigation solution ===========
% It is necessary to have at least three subframes (number 1, 2 and 3) to
% find satellite coordinates. Then receiver position can be found too.
% The function requires all 5 subframes, because the tracking starts at
% arbitrary point. Therefore the first received subframes can be any three
% from the 5.
% One subframe length is 6 seconds, therefore we need at least 30 sec long
% record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
% when tracking has started in a middle of a subframe.


if (settings.msToProcess < 36000) || (sum([trackResults.status] ~= '-') < 4)
    % Show the error message and exit
    disp('Record is too short or too few satellites tracked. Exiting!');
    navSolutions = [];
    eph          = [];
    RAIMresults=[]; unreliableSol= []; TGlobalTest=-1; TlocalTest=-1; GlobalThres=-1; remainingSV=-1; SVexcluded=[];
    return
end


% Initialization RAIM output argument
if settings.RAIM.enableRAIM==false
    RAIMresults=[]; unreliableSol=[];
end
TGlobalTest=-1; TlocalTest=-1; GlobalThres=-1; remainingSV=-1; SVexcluded=[];


%% Find preamble start positions ==========================================
[subFrameStart, activeChnList] = findPreambles(trackResults, settings);

%% Decode ephemerides =====================================================
IONOon=settings.useIonoCorr; % Enables the iono corrections if they are present
[eph, TOW, activeChnList]=decodEph(trackResults,activeChnList,subFrameStart,IONOon);

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    RAIMresults=[]; unreliableSol= []; TGlobalTest=-1; TlocalTest=-1; GlobalThres=-1; remainingSV=-1; SVexcluded=[];
    return
end

%% Initialization =========================================================

% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);

% Save the active channel list. The list contains satellites that are
% tracked and have the required ephemeris data. In the next step the list
% will depend on each satellite's elevation angle, which will change over
% time.
readyChnList = activeChnList;

transmitTime = TOW;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################
endprocessingtime=fix((settings.msToProcess - max(subFrameStart)) /settings.navSolPeriod);
unreliableSol=[];
unrelindex=1; % index for unreliable solution counter
navSolutions.TypeofSol=LSorWLS; % flag to enable W/LS solution. =1 -> LS, =0 -> WLS
%% Initialization of current measurement ==================================
%fprintf('\n');
%fprintf('Time:    \n')
hwb = waitbar(0,'W/LS PVT...');

% Initialization for Protection Levels
PL.HRMS=NaN;
PL.VRMS=NaN;
PL.PRMS=NaN;
PL.HSlopeMax=NaN;
PL.VSlopeMax=NaN;
PL.HPL=NaN;
PL.VPL=NaN;
PL.lambda=NaN;

PL.HRMS2=NaN;
PL.VRMS2=NaN;
PL.PRMS2=NaN;
PL.HSlopeMax2=NaN;
PL.VSlopeMax2=NaN;
PL.HPL2=NaN;
PL.VPL2=NaN;
PL.lambda2=NaN;

% Initialization Ionospheric delay
IONOdelay=zeros(length(readyChnList),endprocessingtime);

% Start processing
for currMeasNr = 1:endprocessingtime
    %fprintf('\b\b\b\b\b\b\b\b\b\b\b')
    %fprintf('Time: %4d',currMeasNr)
    %processing=currMeasNr/fix((settings.msToProcess - max(subFrameStart))/settings.navSolPeriod)*100;
    %    fprintf('\b\b\b\b\b\b\b\b\b\b\b')
    %fprintf('%.2f%% ',processing);
    waitbar(currMeasNr/endprocessingtime, ...
        hwb, ...
        ['PVT with W/LS: Completed ',int2str(currMeasNr/endprocessingtime*100), ...
        ' of ', int2str(100), '%']);
    
    % Exclude satellites, that are belove elevation mask
    if settings.enableLAF==true
        activeChnList=readyChnList;
        %end
        activeChnList = intersect(find(satElev(activeChnList) >= settings.elevationMask), ...
            readyChnList);
    else
        %activeChnList=readyChnList;
        
        activeChnList = intersect(find(satElev >= settings.elevationMask), ...
            readyChnList);
    end
    
    % Test the goodness of SQI for specific PRN at current time
    %     if  1==0 && settings.enableLAF==true && length(activeChnList)>4 %&& settings.RAIM.enableRAIM==true
    %         % Rimuovo i satelliti sottosoglia (hard approach)
    %         %       indexNOTExcludedbySQI=SQIchannels([trackResults(activeChnList).PRN],currMeasNr)>settingsMPDD.SQIthreshold; % Lista PRN salvi
    %         %       fakevariable=not(indexNOTExcludedbySQI);%SQIchannels([trackResults(activeChnList).PRN],currMeasNr)<settingsMPDD.SQIthreshold; % PRN eliminati
    %
    %         % Rimuovo, se possibile, fino al x% dei channel attivi (soft approach)
    %         SVeliminati=ceil(30/100*length(activeChnList));
    %         [kmin,worstIndexList]  = SQImin(SQIchannels([trackResults(activeChnList).PRN],currMeasNr),SVeliminati,settingsMPDD.SQIthreshold); % ne elimino al max SVeliminati
    %         indexNOTExcludedbySQI=not(kmin);
    %         fakevariable=kmin;
    %
    %         PRNexcludedbySQI=[trackResults(activeChnList(fakevariable)).PRN];
    %         %if any(PRNexcludedbySQI)==1, PRNexcludedbySQI=11; end
    %
    %         thresh=settingsMPDD.SQIthreshold;
    %         while sum(indexNOTExcludedbySQI)<4 && thresh>0% se ho pochi SV riduco la soglia
    %             thresh=thresh-0.1
    %             indexNOTExcludedbySQI=SQIchannels([trackResults(activeChnList).PRN],currMeasNr)>thresh; % Lista PRN salvi
    %             fakevariable=not(indexNOTExcludedbySQI); % PRN eliminati
    %             PRNexcludedbySQI=[trackResults(activeChnList(fakevariable)).PRN];
    %             %             indexNOTExcludedbySQI=SQIchannels([trackResults(activeChnList).PRN],currMeasNr)>0 %1 % se ne elimina troppi, non elimino nessuno, lascio la RAIM
    %             %             PRNexcludedbySQI=[];
    %         end
    %
    %         if settings.DOPcontrol==true
    %             NewactiveChnList=activeChnList(indexNOTExcludedbySQI);
    %         else
    %             activeChnList=activeChnList(indexNOTExcludedbySQI);
    %         end
    %
    %     end
    
    % Save list of satellites used for position calculation
    navSolutions.channel.PRN(activeChnList, currMeasNr) = ...
        [trackResults(activeChnList).PRN];
    
    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to position (0,0) in the sky
    % plot.
    navSolutions.channel.el(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    navSolutions.channel.az(:, currMeasNr) = ...
        NaN(settings.numberOfChannels, 1);
    
    %% Find pseudoranges ======================================================
    %            navSolutions.channel.rawP(:, currMeasNr) = calculatePseudoranges(...
    %                    trackResults, ...
    %                    subFrameStart + settings.navSolPeriod * (currMeasNr-1), ...
    %                    activeChnList, settings);
    
    %% BEGIN THUAN
    if currMeasNr>1
        settings.startOffset=settings.startOffset-1000*(navSolutions.dt(currMeasNr-1))/settings.c;
    end
    
    navSolutions.channel.rawP(:, currMeasNr) = calculatePseudoranges(...
        trackResults, ...
        subFrameStart+settings.navSolPeriod*(currMeasNr-1)-1, ...
        activeChnList, settings);
    
    xx=subFrameStart+settings.navSolPeriod*(currMeasNr-1);
    pseudo_correction=zeros(activeChnList(end),1);
    for channelNr = 1:activeChnList(end)
        pseudo_correction(channelNr)=(1023-trackResults(channelNr).remCodePhase(xx(channelNr)-1))/1.023e6*settings.c;
    end
    
    % Correggo le lunghezze dei vettori se non sono uguali
    if length(navSolutions.channel.rawP(:, currMeasNr))==length(pseudo_correction)
        navSolutions.channel.rawP(:, currMeasNr)=navSolutions.channel.rawP(:, currMeasNr)+pseudo_correction;
    else
        difflen=abs(length(navSolutions.channel.rawP(:, currMeasNr))-length(pseudo_correction));
        pseudo_correctionCorr=pseudo_correction;
        for iii=1:difflen
            pseudo_correctionCorr=[pseudo_correctionCorr; Inf];
        end
        navSolutions.channel.rawP(:, currMeasNr)=navSolutions.channel.rawP(:, currMeasNr)+pseudo_correctionCorr;
    end
    %% END  THUAN
    
    
    
    %% Find satellites positions and clocks corrections =======================
    %     [satPositions, satClkCorr] = satpos(transmitTime, ...
    %         [trackResults(activeChnList).PRN], ...
    %         eph, settings);
    
    [satPositions, satClkCorr] = satpos(transmitTime, ...
        [trackResults(:).PRN], ...
        eph, settings);
    
    
    
    %% Find receiver position =================================================
    
    % 3D receiver position can be found only if signals from more than 3
    % satellites are available
    if size(activeChnList,1)>1
        activeChnList=activeChnList';
    end;
    if size(activeChnList,2) > 3
        %% DOP control
        %         if settings.enableLAF==true && settings.DOPcontrol==true && length(NewactiveChnList)<length(activeChnList)
        %             [newChnList,excludedChn]=checkingDOP(satPositions(:,:),navSolutions.channel.rawP(:, currMeasNr), satClkCorr(:),...
        %                 settings,navSolutions.TypeofSol,worstIndexList,activeChnList,NewactiveChnList);
        %             activeChnList=newChnList;
        %             PRNexcludedbySQI=[trackResults(excludedChn).PRN];
        %         end
        
        %%
        Initial_cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty); % .CN0.CNo_SNV(currMeasNr) % navSolutions.channel.PRN(:, currMeasNr),
        cov_noise=Initial_cov_noise;
        
        %         activeChnList1=activeChnList([1:2-1 2+1:end]);
        %         activeChnList=activeChnList1([1:3-1 3+1:end]);
        %         cov_noise=cov_noise(activeChnList,activeChnList);
        %=== Calculate receiver position ==================================
        if settings.enableLAF==true && 1==0
            for SVexclIndex=1:length(PRNexcludedbySQI) %SALVO I PRN ELIMINATI PRIMA DELLA RAIM
                SVexcluded(currMeasNr).PRN(SVexclIndex)=PRNexcludedbySQI(SVexclIndex);
            end
            if isempty(SVexclIndex)==1 % NON HO ELIMINATO SATELLITI
                SVexclIndex=0;
            end
            SVexclIndex=SVexclIndex+1;
            SVexcluded(currMeasNr).PRN(SVexclIndex)=-1; % separo i SV eliminati prima della RAIM con uno -1
            SVexclIndex=SVexclIndex+1;
        else
            SVexclIndex=1;
        end
        
        % W/LS solution
        all_is_ok=-1;
        RecomputePVT=0;
        switcherflag=1; % regola lo switch per SQI exclusion
        repeatPVTflag=true;
        repeatforIONO=0; % parameter for Ionospheric corrections
        
        %%
        while (repeatPVTflag==true)
            while(all_is_ok==-1 && length(activeChnList)>=4)
                
                
                %while repeatforIONO<2 % IONO loop
                    if navSolutions.TypeofSol==1
                        % LS solution
                        [xyzdt, ...
                            navSolutions.channel.el(activeChnList, currMeasNr), ...
                            navSolutions.channel.az(activeChnList, currMeasNr), ...
                            navSolutions.DOP(:, currMeasNr),...
                            ausiliaryVariableforRAIM, IONOdelay(activeChnList,currMeasNr)] = ...
                            leastSquarePos(satPositions(:,activeChnList), ...
                            navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr(activeChnList) * settings.c, ...
                            settings,...
                            eph(trackResults(activeChnList(1)).PRN), transmitTime);
                    else
                        % WLS solution
                        [xyzdt, ...
                            navSolutions.channel.el(activeChnList, currMeasNr), ...
                            navSolutions.channel.az(activeChnList, currMeasNr), ...
                            navSolutions.DOP(:, currMeasNr),navSolutions.WDOP(:, currMeasNr),...
                            ausiliaryVariableforRAIM] = WLS(satPositions(:,activeChnList), ...
                            navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr(activeChnList) * settings.c, ...
                            settings,cov_noise,...
                            eph(trackResults(activeChnList(1)).PRN), transmitTime);
                    end
                    
                    
%                     if eph(trackResults(activeChnList(1)).PRN).iflag==1 % if Ionospheric corrections exist
%                         Alpha=[eph(trackResults(activeChnList(1)).PRN).alpha_0 eph(trackResults(activeChnList(1)).PRN).alpha_1 eph(trackResults(activeChnList(1)).PRN).alpha_2 eph(trackResults(activeChnList(1)).PRN).alpha_3]';
%                         Beta=[eph(trackResults(activeChnList(1)).PRN).beta_0 eph(trackResults(activeChnList(1)).PRN).beta_1 eph(trackResults(activeChnList(1)).PRN).beta_2 eph(trackResults(activeChnList(1)).PRN).beta_3]';
%                         
%                         [IONOcorr]=Ionospheric_Klobuchar_correction(xyzdt(1:3),satPositions(:,activeChnList)',Alpha,Beta,transmitTime,settings.c);
%                         navSolutions.channel.rawP(activeChnList, currMeasNr)=navSolutions.channel.rawP(activeChnList, currMeasNr)+IONOcorr;
%                         repeatforIONO=repeatforIONO+1;
%                     else
%                         repeatforIONO=2; % exit from the IONO loop
%                     end
                    
                %end
                
                RecomputePVT=RecomputePVT+1;
                %%
                if settings.RAIM.enableRAIM==true
                    RAIMresults(currMeasNr)=ausiliaryVariableforRAIM;
                    %cov_noise=PseudoVarianceComputation(navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr(activeChnList) * settings.c,SQI);
                end
                
                if length(activeChnList)>=6 && settings.RAIM.enableRAIM==true && RecomputePVT<settings.RAIM.RecomputePVT %(Con questa condizione abilitata, forzo ad eliminarne solo uno con la RAIM) % mi chiedo se ho abbastanza satelliti per la RAIM
                    % Modifico la lista dei satelliti per ricalcolare la pvt
                    
                    [all_is_ok, measurIndtoexclude, TGlobalTest(currMeasNr),TlocalTest(currMeasNr),GlobalThres,PL]=FaultDetection(xyzdt,RAIMresults(currMeasNr),settings,cov_noise,SQIchannels([trackResults(activeChnList).PRN],currMeasNr));
                    %all_is_ok=1; % Se abilito non faccio Exclusion
                    
                    % Exclusion
                    if all_is_ok<1
                        
                        if all_is_ok==-1 % SV exclusion from the active channel list
                            
                            SVexcluded(currMeasNr).PRN(SVexclIndex)=navSolutions.channel.PRN(activeChnList(measurIndtoexclude),currMeasNr);
                            SVexclIndex=SVexclIndex+1;
                            
                            if measurIndtoexclude==1 % se la misura � da escludere � la prima
                                activeChnList=activeChnList(2:end);
                                
                            else
                                
                                if measurIndtoexclude==length(activeChnList)
                                    activeChnList=activeChnList(1:end-1);
                                else
                                    activeChnList=activeChnList([1:measurIndtoexclude-1 measurIndtoexclude+1:end]);
                                end
                            end
                            
                            %cov_noise=Initial_cov_noise(activeChnList,activeChnList);
                            cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty);
                            %cov_noise=covNoiseMakerEl(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty,navSolutions.channel.el(activeChnList, currMeasNr));

                        else
                            
                            if settings.enableSQI==true % && all_is_ok==0  % SQI could find the fault and avoids solution unreliable
                                [~, minpos]=min(SQIchannels([trackResults(activeChnList).PRN],currMeasNr));
                                if minSQI<settingsMPDD.SQIthreshold
                                    additionalPRNexcluded=trackResults(activeChnList(minpos)).PRN;
                                    
                                    SVexcluded(currMeasNr).PRN(SVexclIndex)=-2;
                                    SVexclIndex=SVexclIndex+1;
                                    SVexcluded(currMeasNr).PRN(SVexclIndex)=additionalPRNexcluded;
                                    SVexclIndex=SVexclIndex+1;
                                    
                                    activeChnList=setdiff(activeChnList,activeChnList(minpos));
                                    RecomputePVT=RecomputePVT-1;
                                    cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty);
                                    all_is_ok=-1;
                                end
                            else
                                unreliableSol(unrelindex)=currMeasNr;
                                unrelindex=unrelindex+1;
                            end
                            
                            
                        end
                    else
                        SVexcluded(currMeasNr).PRN(SVexclIndex)=0;
                    end % if all_is_ok<1
                    
                    %%%%%%%%%%%%%%all_is_ok=1;
                else
                    all_is_ok=1; % se ho meno di 6 satelliti esco dal loop, non faccio la RAIM
                end % length(activeChnList)>=6 && settings.RAIM.enableRAIM==true && RecomputePVT<settings.RAIM.RecomputePVT
            end % end PVT and RAIM
            
            % Possible additional exclusion by SQI:
            % SQI controlled exclusion
            if settings.SQIexclusion==true && length(activeChnList)>=6
                switch switcherflag
                    case 1
                        [minSQI, minpos]=min(SQIchannels([trackResults(activeChnList).PRN],currMeasNr));
                        if minSQI<settingsMPDD.SQIthreshold
                            additionalPRNexcluded=trackResults(activeChnList(minpos)).PRN;
                            
                            SVexcluded(currMeasNr).PRN(SVexclIndex)=-1;
                            SVexclIndex=SVexclIndex+1;
                            SVexcluded(currMeasNr).PRN(SVexclIndex)=additionalPRNexcluded;
                            SVexclIndex=SVexclIndex+1;
                            
                            %Tg=TGlobalTest(currMeasNr);
                            activeChnList_prev=activeChnList;
                            activeChnList=setdiff(activeChnList,activeChnList(minpos));
                            DOP_prev=navSolutions.DOP(:, currMeasNr);
                            switcherflag=switcherflag+1; % becomes "case 2"
                            % re-initialization
                            all_is_ok=-1;
                            RecomputePVT=0;
                            cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty);
                        else
                            repeatPVTflag=false;
                        end
                        
                    case 2
                        %if TGlobalTest(currMeasNr)>Tg %&& Tg>0% reduction of the global test? si per forza, un test delle balle per ora
                        if ((navSolutions.DOP(:, currMeasNr)-DOP_prev)>settings.GDOPbudget*DOP_prev)
                            activeChnList=activeChnList_prev;
                            switcherflag=switcherflag+1;
                            %re-initialization
                            all_is_ok=-1;
                            RecomputePVT=0;
                            cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty);
                            
                        else
                            switcherflag=1;
                            repeatPVTflag=false;
                            %re-initialization
                            all_is_ok=-1;
                            RecomputePVT=0;
                            cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList),1,1,SQIchannels([trackResults(activeChnList).PRN],currMeasNr),1,settings.SQIpenalty);
                            
                        end
                        
                    case 3
                        repeatPVTflag=false;
                        
                end
            else
                repeatPVTflag=false;
            end
        end
        remainingSV(currMeasNr)=length(activeChnList);
        
        %RAIMresults(currMeasNr).chlist=activeChnList;
        %--- Save results -------------------------------------------------
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);
        navSolutions.dt(currMeasNr) = xyzdt(4);
        
        % Update the satellites elevations vector
        satElev(activeChnList) = navSolutions.channel.el(activeChnList, currMeasNr);
        
        
        %=== Correct pseudorange measurements for clocks errors ===========
        navSolutions.channel.correctedP(activeChnList, currMeasNr) = ...
            navSolutions.channel.rawP(activeChnList, currMeasNr) + ...
            satClkCorr(activeChnList)' * settings.c + navSolutions.dt(currMeasNr);
        
        %% Coordinate conversion ==================================================
        
        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr), ...
            navSolutions.height(currMeasNr)] = cart2geo(...
            navSolutions.X(currMeasNr), ...
            navSolutions.Y(currMeasNr), ...
            navSolutions.Z(currMeasNr), ...
            5);
        
        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
            navSolutions.longitude(currMeasNr));
        
        [navSolutions.E(currMeasNr), ...
            navSolutions.N(currMeasNr), ...
            navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
            xyzdt(3), ...
            navSolutions.utmZone);
        
    else % if size(activeChnList, 2) > 3
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
            ': Not enough information for position solution.']);
        
        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;
        
        navSolutions.channel.az(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        navSolutions.channel.el(activeChnList, currMeasNr) = ...
            NaN(1, length(activeChnList));
        
        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rasing
        % satellites will be not included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.
        
    end % if size(activeChnList, 2) > 3
    
    %=== Update the transmit time ("measurement time") ====================
    navSolutions.TOW(currMeasNr)=transmitTime;
    
    transmitTime = transmitTime + settings.navSolPeriod / 1000;
    
    navSolutions.PL(currMeasNr)=PL;
    
end %for currMeasNr...


% Close the waitbar
close(hwb)
fprintf('\n');

end