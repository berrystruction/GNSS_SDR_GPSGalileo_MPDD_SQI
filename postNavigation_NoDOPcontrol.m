function [navSolutions, eph, RAIMresults, unreliableSol,TGlobalTest,TlocalTest,GlobalThres,remainingSV,SVexcluded] = postNavigation(trackResults,settings,settingsMPDD,SQIchannels)
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

for channelNr = activeChnList
    
    %=== Convert tracking output to navigation bits =======================
    
    %--- Copy 5 sub-frames long record from tracking output ---------------
    navBitsSamples = trackResults(channelNr).I_P(subFrameStart(channelNr) - 20 : ...
        subFrameStart(channelNr) + (1500 * 20) -1)';
    
    %--- Group every 20 vales of bits into columns ------------------------
    navBitsSamples = reshape(navBitsSamples, ...
        20, (size(navBitsSamples, 1) / 20));
    
    %--- Sum all samples in the bits to get the best estimate -------------
    navBits = sum(navBitsSamples);
    
    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0);
    
    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    
    %=== Decode ephemerides and TOW of the first sub-frame ================
    [eph(trackResults(channelNr).PRN), TOW] = ...
        ephemeris(navBitsBin(2:1501)', navBitsBin(1));
    
    %--- Exclude satellite if it does not have the necessary nav data -----
    if (isempty(eph(trackResults(channelNr).PRN).IODC) || ...
            isempty(eph(trackResults(channelNr).PRN).IODE_sf2) || ...
            isempty(eph(trackResults(channelNr).PRN).IODE_sf3))
        
        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
    end
end

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
%readyChnList = [1 2 3 4 5];

transmitTime = TOW;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################
endprocessingtime=fix((settings.msToProcess - max(subFrameStart)) /settings.navSolPeriod);
%SVexcluded=[];
unreliableSol=[];
unrelindex=1; % index for unreliable solution counter
%% Initialization of current measurement ==================================
for currMeasNr = 1:endprocessingtime
    
    processing=currMeasNr/fix((settings.msToProcess - max(subFrameStart))/settings.navSolPeriod)*100;
    fprintf('%.2f%% ',processing);
    % Exclude satellites, that are belove elevation mask
    if settings.enableLAF==true
        activeChnList=readyChnList;
        %end
        activeChnList = intersect(find(satElev(activeChnList) >= settings.elevationMask), ...
            readyChnList);
    else
        activeChnList=readyChnList;
        
        activeChnList = intersect(find(satElev >= settings.elevationMask), ...
            readyChnList);
    end
    
    % Test the goodness of SQI for specific PRN at current time
    if settings.enableLAF==true && length(activeChnList)>4%&& settings.RAIM.enableRAIM==true
        % Rimuovo i satelliti sottosoglia (hard approach)
        %       indexNOTExcludedbySQI=SQIchannels([trackResults(activeChnList).PRN],currMeasNr)>settingsMPDD.SQIthreshold; % Lista PRN salvi
        %       fakevariable=not(indexNOTExcludedbySQI);%SQIchannels([trackResults(activeChnList).PRN],currMeasNr)<settingsMPDD.SQIthreshold; % PRN eliminati
        
        % Rimuovo, se possibile, fino al x% dei channel attivi (soft approach)
        SVeliminati=ceil(30/100*length(activeChnList));
        [kmin]  = SQImin(SQIchannels([trackResults(activeChnList).PRN],currMeasNr),SVeliminati,settingsMPDD.SQIthreshold); % ne elimino al max SVeliminati
        indexNOTExcludedbySQI=not(kmin);
        fakevariable=kmin;
        
        PRNexcludedbySQI=[trackResults(activeChnList(fakevariable)).PRN];
        %if any(PRNexcludedbySQI)==1, PRNexcludedbySQI=11; end
        
        thresh=settingsMPDD.SQIthreshold;
        while sum(indexNOTExcludedbySQI)<4 && thresh>0% se ho pochi SV riduco la soglia
            thresh=thresh-0.1
            indexNOTExcludedbySQI=SQIchannels([trackResults(activeChnList).PRN],currMeasNr)>thresh; % Lista PRN salvi
            fakevariable=not(indexNOTExcludedbySQI); % PRN eliminati
            PRNexcludedbySQI=[trackResults(activeChnList(fakevariable)).PRN];
            %             indexNOTExcludedbySQI=SQIchannels([trackResults(activeChnList).PRN],currMeasNr)>0 %1 % se ne elimina troppi, non elimino nessuno, lascio la RAIM
            %             PRNexcludedbySQI=[];
        end
        activeChnList=activeChnList(indexNOTExcludedbySQI);
    end
    
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
    %        navSolutions.channel.rawP(:, currMeasNr) = calculatePseudoranges(...
    %                trackResults, ...
    %                subFrameStart + settings.navSolPeriod * (currMeasNr-1), ...
    %                activeChnList, settings);
    
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
        
        
        Initial_cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList)); % .CN0.CNo_SNV(currMeasNr) % navSolutions.channel.PRN(:, currMeasNr),
        cov_noise=Initial_cov_noise;
        
        %         activeChnList1=activeChnList([1:2-1 2+1:end]);
        %         activeChnList=activeChnList1([1:3-1 3+1:end]);
        %         cov_noise=cov_noise(activeChnList,activeChnList);
        %=== Calculate receiver position ==================================
        if settings.enableLAF==true
            for SVexclIndex=1:length(PRNexcludedbySQI) %SALVO I PRN ELIMINATI PRIMA DELLA RAIM
                SVexcluded(currMeasNr).PRN(SVexclIndex)=PRNexcludedbySQI(SVexclIndex);
            end
            SVexclIndex=SVexclIndex+1;
            SVexcluded(currMeasNr).PRN(SVexclIndex)=-1; % separo i SV eliminati prima della RAIM con uno -1
            SVexclIndex=SVexclIndex+1;
        else
            SVexclIndex=1;
        end
        % LS solution
        all_is_ok=-1;
        
        while(all_is_ok==-1 && length(activeChnList)>=4)
            % LS solution
            [xyzdt, ...
                navSolutions.channel.el(activeChnList, currMeasNr), ...
                navSolutions.channel.az(activeChnList, currMeasNr), ...
                navSolutions.DOP(:, currMeasNr),...
                ausiliaryVariableforRAIM] = ...
                leastSquarePos(satPositions(:,activeChnList), ...
                navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr(activeChnList) * settings.c, ...
                settings);
            
            if settings.RAIM.enableRAIM==true
                RAIMresults(currMeasNr)=ausiliaryVariableforRAIM;
            end
            % WLS solution
            %             [xyzdt, ...
            %                 navSolutions.channel.el(activeChnList, currMeasNr), ...
            %                 navSolutions.channel.az(activeChnList, currMeasNr), ...
            %                 navSolutions.DOP(:, currMeasNr),...
            %                 RAIMresults(currMeasNr)] = WLS(satPositions(:,activeChnList), ...
            %                 navSolutions.channel.rawP(activeChnList, currMeasNr)' + satClkCorr(activeChnList) * settings.c, ...
            %                 settings,cov_noise);
            
            if length(activeChnList)>6 && settings.RAIM.enableRAIM==true% mi chiedo se ho abbastanza satelliti per la RAIM
                % Modifico la lista dei satelliti per ricalcolare la pvt
                [all_is_ok, measurIndtoexclude, TGlobalTest(currMeasNr),TlocalTest(currMeasNr),GlobalThres]=FaultDetection(RAIMresults(currMeasNr),settings,cov_noise);
                % Exclusion
                if all_is_ok<1
                    
                    if all_is_ok==-1 % SV exclusion from the active channel list
                        
                        SVexcluded(currMeasNr).PRN(SVexclIndex)=navSolutions.channel.PRN(activeChnList(measurIndtoexclude),currMeasNr);
                        SVexclIndex=SVexclIndex+1;
                        
                        if measurIndtoexclude==1 % se la misura è da escludere è la prima
                            activeChnList=activeChnList(2:end);
                            
                        else
                            
                            if measurIndtoexclude==length(activeChnList)
                                activeChnList=activeChnList(1:end-1);
                            else
                                activeChnList=activeChnList([1:measurIndtoexclude-1 measurIndtoexclude+1:end]);
                            end
                        end
                        
                        %cov_noise=Initial_cov_noise(activeChnList,activeChnList);
                        cov_noise=covNoiseMaker(trackResults(activeChnList),currMeasNr,length(activeChnList));
                    else
                        unreliableSol(unrelindex)=currMeasNr;
                        unrelindex=unrelindex+1;
                    end
                else
                    SVexcluded(currMeasNr).PRN(SVexclIndex)=0;
                end % if all_is_ok<1
                
                all_is_ok=1; % Con questa condizione abilitata, forzo ad eliminarne solo uno con la RAIM
            else
                all_is_ok=1; % se ho meno di 6 satelliti esco dal loop, non faccio la RAIM
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
        %        if settings.enableLAF==true
        satElev(activeChnList) = navSolutions.channel.el(activeChnList, currMeasNr);
        %         else
        %             satElev= navSolutions.channel.el(:, currMeasNr);
        %         end
        
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
    transmitTime = transmitTime + settings.navSolPeriod / 1000;
    
end %for currMeasNr...


fprintf('\n');

end