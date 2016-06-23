function [navSolutions,RAIMresults,unreliableSol,TGlobalTest,TlocalTest,GlobalThres,remainingSV,SVexcluded,Tk,ThrTk,PVTrunmode]=PVT_Kalman_FallBack(trackResults,settings,kalmansettings,LSorWLS,SQI,SQIthreshold)
%% FILE to LOAD
%mat_file = 'trackingResultsL1_morning26.mat';
%load(mat_file);

channels_GPS_L1 = trackResults;
%clear trackResults;
%clear acqResults;
navSolutions.TypeofSol=LSorWLS;
%SimTime =settings.msToProcess/1000;%10*60; % in s


%% ATTENTION!!! This overwrites the content of setting in TrackResults!!
settings.enable_smoothL1    = 0;
settings.smoothL1_weights   = 100;
settings.enableLSfallback   = 1;
settings.enableSQI          = 0;

settingsL1 = settings;

%% Initialization
startKalman = kalmansettings.startKalman;
f_sampling = settings.samplingFreq; % sampling frequency [Hz]
sol = settings.c;

%if settings.RAIM.enableRAIM==false
%    RAIMresults=[]; unreliableSol=[];
%end
TGlobalTest=-1; TlocalTest=-1; GlobalThres=-1; remainingSV=-1; SVexcluded=[];Tk=-1;

% PVT rate
%settings.navSolPeriod = 1000; % in ms
Rate = settings.navSolPeriod/1000; % in seconds
num_samples_Rate = Rate*settings.samplingFreq; % in samples
ConstantsInit;

% *********************************************************************************
% *                          MAIN                                                 *
% *********************************************************************************

% list of channels (not PRN id) involved in PVT computation
% Note:
%       1) For single system, the number of channels of that system must be
%       larger than 4. For combined solution, both systems are used, the
%       number of channels from the two systems must be larger than 5.

%%------ ELIMINATION CHANNELS BACATS -----------------------------------
%      channels_GPS_L1(5) = [];
%      channels_GPS_L1(1) = [];

activeChnList_GPS_L1 = 1:length(find([channels_GPS_L1(:).status]=='T'));

%% Data demodulation
[channels_GPS_L1,subFrameStart] = GPS_data_dem(channels_GPS_L1,settings.T_int,settings); % ORIGINAL: settings.L1Tc - "upsample"?????

for channelNr = activeChnList_GPS_L1
    if (isempty(channels_GPS_L1(channelNr).eph.IODC) || ...
            isempty(channels_GPS_L1(channelNr).eph.IODE_sf2) || ...
            isempty(channels_GPS_L1(channelNr).eph.IODE_sf3))
        
        %--- Exclude channel from the list (from further processing) ------
        activeChnList_GPS_L1 = setdiff(activeChnList_GPS_L1, channelNr);
    end
end

% After data demodulation the structure of the channels increase with the following variables:
%  channels[index].WordStart: location in symbol of the starting of the
%                            subframe for GPS;
%                 .eph: ephemeris demodulated from the data symbol

% location measured in symbol indicates where the subframe start
[startClock_PVT_GPS_L1,channels_GPS_L1] = startPVT_GPS(channels_GPS_L1, num_samples_Rate, settingsL1);
KGPS = initKalmanGPS;

% For calculating pseudorange
startOffset_d0_GPS = settings.startOffset;
startOffset_d1_GPS = settings.startOffset;
startOffset_d2_GPS = settings.startOffset;
startOffset_d0_GPS_L1 = settings.startOffset;

%length_time = SimTime/Rate;
length_time=fix((settings.msToProcess - max(subFrameStart)) /settings.navSolPeriod);


N_sat = length(activeChnList_GPS_L1);
bias  = 0; % for carrier smoothing
pseudo_PhaseL1     =  zeros(N_sat,length_time);
pseudo_smoothedL1  =  zeros(N_sat,length_time);
pseudo_codeL1      =  zeros(N_sat,length_time);

for index = 1:numel(channels_GPS_L1)
    channels_GPS_L1(index).ind = 1;
end


% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);


%% MAIN LOOP ---------------------------------------------------------------------
loop = 'true';
fprintf('\n')
%fprintf('Time:    \n')
pos(8) = 0;

unrelindex=1; % index for unreliable solution counter
unreliableSol=[];

% Initialization for Protection Levels
PL.HRMS=NaN;
PL.VRMS=NaN;
PL.PRMS=NaN;
PL.HSlopeMax=NaN;
PL.VSlopeMax=NaN;
PL.HPL=NaN;
PL.VPL=NaN;
PL.lambda=NaN;

% Init for LS fallback

LSfallback=0; %flag for LS fallback. Initially set to 0 to disable, will turn on in case of problem
LSfallbackcounter=1; %a counter of epochs that run LS fallback
LSfallbackround=5; %maximum epochs that run by LS fallback
thresholdcoef=2.7;
PVTrunmode=zeros(1,length_time);
avgwindows=8;


readyChnList = activeChnList_GPS_L1;
hwb = waitbar(0,'Kalman PVT...');

for currMeasNr = 1:length_time
    
    if currMeasNr < startKalman
        cc= ['PVT with Least Squares: Completed ',int2str(currMeasNr/length_time*100), ' of ', int2str(100), '%'];
    else
        cc= ['PVT with Kalman Filter: Completed ',int2str(currMeasNr/length_time*100), ' of ', int2str(100), '%'];
    end
    
    waitbar(currMeasNr/length_time,hwb,cc);
    
    % Exclude satellites, that are belove elevation mask
    if settings.enableLAF==true
        activeChnList_GPS_L1=readyChnList;
        %end
        activeChnList_GPS_L1 = intersect(find(satElev(activeChnList_GPS_L1) >= settings.elevationMask), ...
            readyChnList);
    else
        activeChnList_GPS_L1=readyChnList;
        
        activeChnList_GPS_L1 = intersect(find(satElev >= settings.elevationMask), ...
            readyChnList);
    end
    
    
    
    %% Location Computation
    Clock = startClock_PVT_GPS_L1 + num_samples_Rate*(currMeasNr-1);  % the clock_time where to compute the PVT is increased according to the PVT rate
    [channels_GPS_L1,Doppler_GPS_L1] = Counter_clock_Update_GPS(channels_GPS_L1, Clock, settingsL1, currMeasNr, pos(8));
    
    [transmitTime_increment_GPS_L1] = Compute_SatTX_time_GPS(channels_GPS_L1, settingsL1);
    if currMeasNr > 1
        startOffset_d0_GPS = (RX_GPS_tow(1) - max(transmitTime_increment_GPS_L1)) * 1000; % in ms
    end
    
    [Timediff_GPS_L1, TravelTime_GPS_L1, raw_pseudoranges_GPS_L1, delta_pseudoranges_GPS_L1] = ...
        sample2time_GPS(channels_GPS_L1, Doppler_GPS_L1, startOffset_d0_GPS, settingsL1);
    
    %% Satellite position
    [channels_GPS_L1] = GPS_satpos(channels_GPS_L1, activeChnList_GPS_L1, transmitTime_increment_GPS_L1);
    if (Clock > channels_GPS_L1(1,1).absoluteSample(end)) %check if we have reached the end of the file
        loop = 'false';
        break;
    end
    pseudoranges1_GPS_L1 = raw_pseudoranges_GPS_L1(activeChnList_GPS_L1)' + ...
        [channels_GPS_L1(activeChnList_GPS_L1).satClkCorr]*sol; % corrected pseudorange with clock errors
    
    i1 = 1;
    satpos = zeros(3,length(activeChnList_GPS_L1));
    satvel = satpos;
    for i11 = activeChnList_GPS_L1
        satpos(:,i1) = channels_GPS_L1(i11).satPositions;
        satvel(:,i1) = channels_GPS_L1(i11).satVelocity;
        i1 = i1+1;
    end
    
    if (settings.enable_smoothL1)
        %% CARRIER SMOOTHING
        for i11 = activeChnList_GPS_L1
            channels_GPS_L1(i11).Lpc(1,currMeasNr) = channels_GPS_L1(i11).Lpc(1,currMeasNr) + bias; % tolgo il bias anche sulle misure di fase
            pseudo_codeL1 (i11,currMeasNr) = pseudoranges1_GPS_L1(i11);
            
            if (currMeasNr==2)
                pseudo_PhaseL1(i11,currMeasNr) = pseudo_codeL1 (i11,currMeasNr-1) - (channels_GPS_L1(i11).Lpc(currMeasNr)...
                    - channels_GPS_L1(i11).Lpc(currMeasNr-1));
            else
                if (currMeasNr~=1)
                    pseudo_PhaseL1(i11,currMeasNr) = pseudo_smoothedL1 (i11,currMeasNr-1) - (channels_GPS_L1(i11).Lpc(currMeasNr)...
                        - channels_GPS_L1(i11).Lpc(currMeasNr-1));
                end
            end
            pseudo_smoothedL1(i11,currMeasNr) = 1/settings.smoothL1_weights * pseudo_codeL1(i11,currMeasNr) + ...
                (settings.smoothL1_weights-1) / settings.smoothL1_weights * pseudo_PhaseL1(i11,currMeasNr);
        end
    end
    
    SVexclIndex=1; % index of excluded satellites
    % Save list of satellites used for position calculation
    navSolutions.channel.PRN(activeChnList_GPS_L1, currMeasNr) =[channels_GPS_L1(activeChnList_GPS_L1).PRN];
    
    if currMeasNr < startKalman || LSfallback>0 % || (currMeasNr > 106 && currMeasNr < 135) %forcefully fallback to LS
        % compute the user's position using LS
        if(settings.enable_smoothL1 == 1)
            if(currMeasNr == 1)
                [pos,satElev(activeChnList_GPS_L1),az,dop] = leastSquarePos_Doppler(satpos,pseudo_codeL1(:,currMeasNr),channels_GPS_L1(activeChnList_GPS_L1), Doppler_GPS_L1, settings);
            else
                % [pos,satElev(activeChnList_GPS),az,dop] = leastSquarePos_flex(satpos,pseudo_sm(:,currMeasNr),channels_GPS(activeChnList_GPS),Doppler_GPS,transmitTime_increment_GPS,pseudo_flag);
                [pos,satElev(activeChnList_GPS_L1),az,dop] = leastSquarePos_Doppler(satpos,pseudo_smoothedL1(:,currMeasNr),channels_GPS_L1(activeChnList_GPS_L1), Doppler_GPS_L1, settings);
            end
            bias = bias+pos(4)-pos(8);
            
        else % No carrier smoothing
            
            all_is_ok=-1;
            RecomputePVT=0;
            %cov_noise  = covNoiseMaker(channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,length(activeChnList_GPS_L1));
            while(all_is_ok==-1 && length(activeChnList_GPS_L1)>=4)
                cov_noise  = covNoiseMaker(channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,length(activeChnList_GPS_L1),...
                    pseudoranges1_GPS_L1(activeChnList_GPS_L1),0,SQI([channels_GPS_L1(activeChnList_GPS_L1).PRN],currMeasNr),SQIthreshold,settings.SQIpenalty);
                % Navigation solution computation
                if LSorWLS==1
                    % LS
                    [pos, satElev(activeChnList_GPS_L1),az,dop,RAIMresult] =...
                        leastSquarePos_Doppler(satpos(:,activeChnList_GPS_L1), pseudoranges1_GPS_L1(activeChnList_GPS_L1)', channels_GPS_L1(activeChnList_GPS_L1), Doppler_GPS_L1(activeChnList_GPS_L1), settings);
                else
                    % WLS
                    [pos, satElev(activeChnList_GPS_L1),az,dop,RAIMresult,wdop] =...
                        weightedLeastSquarePos_Doppler(satpos(:,activeChnList_GPS_L1), pseudoranges1_GPS_L1(activeChnList_GPS_L1)', channels_GPS_L1(activeChnList_GPS_L1), Doppler_GPS_L1(activeChnList_GPS_L1), settings,cov_noise);

                    navSolutions.WDOP(:, currMeasNr)       = wdop;
                end
                
                
                RecomputePVT=RecomputePVT+1;
                RAIMresults(currMeasNr) = RAIMresult;
                
                % ... RAIM
                if length(activeChnList_GPS_L1)>=6 && settings.RAIM.enableRAIM==true && RecomputePVT<3
                    %                     RAIMresults(currMeasNr) = RAIMresult;
                    [all_is_ok, measurIndtoexclude, TGlobalTest(currMeasNr,1),TlocalTest(currMeasNr),GlobalThres,PL...Tk(currMeasNr),ThrTk(currMeasNr)
                        ]=FaultDetection(pos,RAIMresults(currMeasNr),settings,cov_noise);
                    %all_is_ok=1; % To avoid the exclusion
                    % Exclusion
                    ExclusionSatellites
                    if currMeasNr==1
                        M=TGlobalTest(currMeasNr);
                        M_old=TGlobalTest(currMeasNr);
                    else
                        M_old=M;
                        M=M_old+(TGlobalTest(currMeasNr)-M_old)/currMeasNr;
                    end
                else
                    all_is_ok=1;% se ho meno di 6 satelliti esco dal loop, non faccio la RAIM
                end % length(activeChnList_GPS_L1)>=6 && settings.RAIM.enableRAIM==true && RecomputePVT<2
                
                
            end
            
        end
        
        navSolutions.channel.el(activeChnList_GPS_L1,currMeasNr)=satElev(activeChnList_GPS_L1);
        navSolutions.channel.az(activeChnList_GPS_L1,currMeasNr)=az;
        navSolutions.channel.PRN(activeChnList_GPS_L1,currMeasNr)=[channels_GPS_L1(activeChnList_GPS_L1).PRN];
        %         navSolutions.DOP(:,currMeasNr)=dop;
        navSolutions.TypeofSol=1; % LS
        posprec=pos(1:3);
        
        LSfallbackcounter=LSfallbackcounter+1;
        if LSfallbackcounter>LSfallbackround %check if fallback has run enough, then switch to Kalman again
            LSfallback=0;
        end
        
        %Update the log
        PVTrunmode(currMeasNr)=1;
        TGlobalTest(currMeasNr,2)=NaN;
        
    else % Kalman Filter branch
        
        
        % [pseudoranges_GPS,el] = pseudo2true_V2(pseudoranges1_GPS', zeros(size(activeChnList_GPS)), pos(1:3), 0, satpos,channels_GPS(activeChnList_GPS)); % from pseudoranges to ranges
        pos(1:3) = pos(1:3) + pos(5:7)*Rate;
        pos(4) = pos(4)-pos(8)*Rate;
        
        % EKF solution
        all_is_ok=-1;
        RecomputePVT=0;
        while(all_is_ok==-1 && length(activeChnList_GPS_L1)>=4)
            
            pseudo=0; % da finire
            KGPS = initSigmaKalmanGPS(KGPS,channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,kalmansettings.staticposition,pseudo,SQI([channels_GPS_L1(activeChnList_GPS_L1).PRN],currMeasNr),SQIthreshold,0);
        
            
            % compute the user's position using KF
            AAA = pseudoranges1_GPS_L1(activeChnList_GPS_L1)/sol * Const.WEDOT; % earth rotation compensation
            
            
            [pseudoranges_GPS_L1,el,az] = pseudo2true(pseudoranges1_GPS_L1(activeChnList_GPS_L1)', zeros(size(activeChnList_GPS_L1)), pos(1:3), 0, satpos(:,activeChnList_GPS_L1),sol); % from pseudoranges to ranges
            %pos(1:3) = pos(1:3) + pos(5:7)*Rate;
            
            % EKF solutions and...
            if kalmansettings.staticposition==false
                [KGPS.stateAPosteriori, KGPS.PAPostCovariance,dop,RAIMparam] = KfilterGPS(KGPS, pseudoranges_GPS_L1, delta_pseudoranges_GPS_L1(activeChnList_GPS_L1), pos(1:3), pos(5:7), satpos(:,activeChnList_GPS_L1), satvel(:,activeChnList_GPS_L1), pos(4), pos(8), AAA);
                pos(1:7) = pos(1:7) + KGPS.stateAPosteriori(1:7);
                pos(8) = pos(8) - KGPS.stateAPosteriori(8); %   pos(end) = pos(end) - 2*KGPS.stateAPosteriori(end);
            
                
%                 for IndexU = 1 : length(pseudoranges_GPS_L1)
%                 %% Compute the Range with the relativistic correction and the satellite clock bias
%                         dX 		= RAIMparam.Hgeometric(IndexU,1) - pos(1);
%                         dY 		= RAIMparam.Hgeometric(IndexU,2) - pos(2);
%                         dZ 		= RAIMparam.Hgeometric(IndexU,3) - pos(3);
% 
%                         Range(IndexU) = sqrt(power(dX,2)+power(dY,2)+power(dZ,2));
% 
%                     % RAIMparam.w(IndexU,1)=(pseudoranges_GPS_L1(IndexU))-(Range(IndexU)+pos(4));
%                      
%                 end
                RAIMparam.w=RAIMparam.y-RAIMparam.H*KGPS.stateAPosteriori(1:4);
                
            else
                [KGPS.stateAPosteriori, KGPS.PAPostCovariance,dop,RAIMparam] = KfilterGPSstatic(KGPS, pseudoranges_GPS_L1, delta_pseudoranges_GPS_L1(activeChnList_GPS_L1), pos(1:3), pos(5:7), satpos(:,activeChnList_GPS_L1), satvel(:,activeChnList_GPS_L1), pos(4), pos(8), AAA,posprec);
                pos = pos + KGPS.stateAPosteriori;
                pos(end) = pos(end) - KGPS.stateAPosteriori(end); %   pos(end) = pos(end) - 2*KGPS.stateAPosteriori(end);
                posprec=pos(1:3);
            end
            KGPS.stateAPosteriori1 = KGPS.stateAPosteriori;
            KGPS.stateAPosteriori = zeros(size(KGPS.stateAPosteriori));
            %KGPS.KFOutput=KFOutput;
            
            RecomputePVT=RecomputePVT+1;
            % ... RAIM
            if length(activeChnList_GPS_L1)>=6 && settings.RAIM.enableRAIM==true && RecomputePVT<3
                RAIMresults(currMeasNr) = RAIMparam;
                [all_is_ok, measurIndtoexclude, TGlobalTest(currMeasNr,2),TlocalTest(currMeasNr),GlobalThres,PL...,Tk(currMeasNr),ThrTk(currMeasNr)
                    ]=FaultDetection(pos,RAIMresults(currMeasNr),settings,...
                    KGPS.RObservationNoiseCovariance(1:length(activeChnList_GPS_L1),1:length(activeChnList_GPS_L1)),KGPS);
                %all_is_ok=1; % To avoid the exclusion
                % Exclusion
                ExclusionSatellites
                
                %---Fallback analysis
                if currMeasNr>avgwindows+1
                    recenttest=TGlobalTest(currMeasNr-avgwindows-1:currMeasNr-1);
                else
                    recenttest=TGlobalTest(1:currMeasNr-1);
                end
                avgtest=mean(recenttest);
                %fallbackthreshold=avgtest*thresholdcoef;
                fallbackthreshold=5;
                
                %fallbackthreshold=M*thresholdcoef;
                if TGlobalTest(currMeasNr,2)>fallbackthreshold && settings.enableLSfallback>0
                    LSfallback=1;
                    LSfallbackcounter=1;
                    %LSfallbackround=LSfallbackround+4;
                    disp('Fallback to LS now!');
                end
                
                %Update the running average
                if currMeasNr==1
                    M=TGlobalTest(currMeasNr);
                    M_old=TGlobalTest(currMeasNr);
                else
                    M_old=M;
                    M=M_old+(TGlobalTest(currMeasNr)-M_old)/currMeasNr;
                end
                
            else
                all_is_ok=1;
            end
            
        end
        
        navSolutions.channel.el(activeChnList_GPS_L1,currMeasNr)=el/pi*180;
        navSolutions.channel.az(activeChnList_GPS_L1,currMeasNr)=az/pi*180;
        navSolutions.channel.PRN(activeChnList_GPS_L1,currMeasNr)=[channels_GPS_L1(activeChnList_GPS_L1).PRN];
        %         end
        
        %Update the log
        PVTrunmode(currMeasNr)=2;
        TGlobalTest(currMeasNr,1)=NaN;
    end
    
    remainingSV(currMeasNr)=length(activeChnList_GPS_L1);
    
    Rx_referenceTOW = transmitTime_increment_GPS_L1 + (raw_pseudoranges_GPS_L1/sol-(pos(4)/sol)); % in seconds
    delta_time = num_samples_Rate/f_sampling; % seconds
    RX_GPS_tow = Rx_referenceTOW + delta_time/(1.0-pos(8)/sol);
    % f_sampling = f_sampling/(1.0-pos(8)/sol);
    % settings_f_sampling = f_sampling;
    
    % Save results
    pseudorangeL1.PRN(activeChnList_GPS_L1, currMeasNr)           = [channels_GPS_L1(activeChnList_GPS_L1).PRN];
    pseudorangeL1.rawP(activeChnList_GPS_L1, currMeasNr)          = raw_pseudoranges_GPS_L1(activeChnList_GPS_L1);
    pseudorangeL1.t(activeChnList_GPS_L1, currMeasNr)             = Clock;
    pseudorangeL1.correctedP(activeChnList_GPS_L1, currMeasNr)    = pseudoranges1_GPS_L1(activeChnList_GPS_L1);
    pseudorangeL1.bias(activeChnList_GPS_L1, currMeasNr)          = pos(4);
    pseudorangeL1.GPStime(activeChnList_GPS_L1, currMeasNr)       = RX_GPS_tow(1);%Rx_referenceTOW(1);
    pseudorangeL1.Doppler_meter(activeChnList_GPS_L1, currMeasNr) = delta_pseudoranges_GPS_L1(activeChnList_GPS_L1);
    
    if(settings.enable_smoothL1 == 1)
        pseudorangeL1.smoothedP(activeChnList_GPS_L1, currMeasNr) = pseudo_smoothedL1 (:,currMeasNr)';
    end
    
    %% Data storage
    ECEF = [pos(1) pos(2) pos(3)];
    [LLH] = ECEF2LLH(ECEF);
    Separation  = WGS84Separation(LLH(1),LLH(2));
    LLH(3)	= LLH(3) - Separation;
    
    navSolutions.X(currMeasNr)             = pos(1);
    navSolutions.Y(currMeasNr)             = pos(2);
    navSolutions.Z(currMeasNr)             = pos(3);
    navSolutions.latitude(currMeasNr)      = LLH(1);
    navSolutions.longitude(currMeasNr)     = LLH(2);
    navSolutions.height(currMeasNr)        = LLH(3);
    navSolutions.Vx(currMeasNr)            = pos(5);
    navSolutions.Vy(currMeasNr)            = pos(6);
    navSolutions.Vz(currMeasNr)            = pos(7);
    navSolutions.V(currMeasNr)              = sqrt(pos(5)^2+pos(6)^2+pos(7)^2);
    navSolutions.Clock_GPS(currMeasNr)     = pos(4);
    navSolutions.Drifr_Clk_GPS(currMeasNr) = pos(8);
    navSolutions.DOP(:,currMeasNr)         = dop;
    navSolutions.channel.rawP(activeChnList_GPS_L1,currMeasNr)=  raw_pseudoranges_GPS_L1(activeChnList_GPS_L1);
    % TOW
    navSolutions.TOW(currMeasNr)= Rx_referenceTOW(1);
    navSolutions.PL(currMeasNr)=PL;
    
    
    % Convert to UTM coordinate system
    navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
        navSolutions.longitude(currMeasNr));
    
    [navSolutions.E(currMeasNr), ...
        navSolutions.N(currMeasNr), ...
        navSolutions.U(currMeasNr)] = cart2utm(pos(1), pos(2), pos(3), ...
        navSolutions.utmZone);
    
end

%%
% Close the waitbar
close(hwb)
%fprintf('\n');
for n=1:length(activeChnList_GPS_L1)
    T_GD_L1(n) = channels_GPS_L1(n).eph.T_GD;
end

%save('pseudorangeL1','pseudorangeL1');

% Write points to a KML file.
lat     = [navSolutions.latitude];
long    = [navSolutions.longitude];
height  = navSolutions.height;
PRN_GPS_L1 = [];
for n=1:length(activeChnList_GPS_L1)
    PRN_GPS_L1 = [PRN_GPS_L1 channels_GPS_L1([activeChnList_GPS_L1(n)]).PRN];
end
%filename2 = 'Nav_Sol_GPS_L1.mat';
%filename = 'Nav_Sol_GPS_L1.kml';

%data_storage(filename2,lat(1:end-1),long(1:end-1),PRN_GPS_L1);

%save('navSolutions.mat','navSolutions');
%KML_googleEarth(myfolder, lat(1:end), long(1:end),height(1:end));

% plotNavigation(navSolutions,settings,myfolder)

if settings.enableLSfallback
    figure;
    plot(PVTrunmode,'o-'), grid on
end 

end
