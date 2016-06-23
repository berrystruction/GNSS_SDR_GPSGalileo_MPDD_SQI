function [trackResults, channel]= tracking(fid,channel,settings,settingsMPDD,FLL_results)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
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
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

codePeriods = (settings.msToProcess);     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, ...
    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
    settings.pllDampingRatio, ...
    0.25); % original value: 0.25);

%% C/No estimator
% initialization
for channelNr = 1:settings.numberOfChannels
    trackResults(channelNr).CN0.CNo_Emanuela = zeros(1,ceil(codePeriods/settingsMPDD.CNo_WINDOW));
    trackResults(channelNr).CN0.CNo_bluebook = zeros(1,ceil(codePeriods/settingsMPDD.CNo_WINDOW));
    trackResults(channelNr).CN0.CNo_SNV = zeros(1,ceil(codePeriods/settingsMPDD.CNo_WINDOW));
    trackResults(channelNr).CN0.CNo_MM = zeros(1,ceil(codePeriods/settingsMPDD.CNo_WINDOW));
    trackResults(channelNr).CN0.CNo_Bea = zeros(1,ceil(codePeriods/settingsMPDD.CNo_WINDOW));
    
    trackResults(channelNr).rem_fraction = zeros(1, settings.msToProcess);
    trackResults(channelNr).rem_integer  = zeros(1, settings.msToProcess);
end

FIFO_IP = zeros(settings.numberOfChannels,settingsMPDD.CNo_WINDOW/(settings.T_int*1e3));
FIFO_QP = zeros(settings.numberOfChannels,settingsMPDD.CNo_WINDOW/(settings.T_int*1e3));

% ausiliary variable
CN0aus.CNo_Emanuela=0;
CN0aus.CNo_bluebook=0;
CN0aus.CNo_SNV=0;
CN0aus.CNo_MM=0;
CN0aus.CNo_Bea=0;



%% LAF initialization
if settings.enableLAF==true
    
    % skp_msec=settings.seek_sec;
    %msec=settings.track_time; msecthreshold=settings.mvgAvgtime;
    M=settingsMPDD.M; Npoints=settingsMPDD.multicorr.Npoints; maxlag=settingsMPDD.multicorr.maxlag;
    version=settingsMPDD.LAFversion;
    
    alignment=true;
    if alignment==true
        Nallignment=2*Npoints+1-floor(Npoints/2);
        %Nallignment=2*Npoints+1;
        N=2*Npoints+1;
    else
        Nallignment=2*Npoints+1;
        N=Nallignment;
    end
    
    %version=2; % if version==2, LAF standard function (covariance), version~=2
    %autocorrelation method of LAF or others
    %if version==2
    %    correctNumber=Nallignment-M+1;% corretta dimensione POST LAF
    %else
    %    correctNumber=Nallignment+M-1;% N+M-1;  autocorrelation method
    %    %correctNumber=Nallignment-M+1+floor(tott); % CLS-LAF
    %end
    
    u1=zeros(N,1); %u1(ceil(N/2))=1;
    d1=zeros(N,1); %d1(ceil(N/2))=1;
    
    stepp=maxlag/Npoints; % parametri per il multicorrelatore
    minlag=stepp;
    
    %%%correctDimension=floor(codePeriods/settingsMPDD.mvgAvgtime/(settings.T_int*1e3)); % dimensione della matrice (temporalmente) corretta
    
    %     powersFinal=zeros(correctDimension,5);  % Powers matrix
    %     yFinal=zeros(correctNumber,correctDimension);
    %     ynewFinal=zeros(correctNumber,correctDimension);
    %     wFinal=zeros(M,correctDimension);         % weights of ideal correlation
    %     error_vecFinal=zeros(M,correctDimension);   % error vector for covariance calculation
    %     eFinal=zeros(correctNumber,correctDimension);
    %     Ufinal=zeros(correctNumber,M,correctDimension);           % Correlation Matrix
    %     Dfinal=zeros(correctNumber,correctDimension);
    %     ufinal=zeros(Nallignment,correctDimension);
    %     dfinal=ufinal;
    %
    %     mvgAvgindex=1; %index for mvgAvg matrices
    
    %%%tau_m=zeros(Nallignment,correctDimension);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hwb = waitbar(0,'Tracking...');
set(hwb,'Name',['Data file ' num2str(settings.datafileIndex) ' of ' num2str(settings.N_of_datafile)]);
%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels
    
    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;
        
        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample)
        
        %         fseek(fid, ...
        %             settings.skipNumberOfBytes + channel(channelNr).codePhase-1, ...
        %             'bof');
        switch settings.GNSS_signal
            case ('GAL_E1b')
                fseek(fid, ...
                    FLL_results.phaseFLL((trackResults(channelNr).PRN))-1, ...
                    'bof');          
            case ('GPS_L1')              
                fseek(fid,(settings.NSample*settings.BytePerSample*(FLL_results.phaseFLL((trackResults(channelNr).PRN))-1)+...
                    round((20-FLL_results.cnt_skp((trackResults(channelNr).PRN)))*settings.samplingFreq*settings.NSample*settings.BytePerSample*1e-3)),'bof');
                % IMPORTANT: mycodephase: number of samples skipped during the acquisition  and processed by the FLL from the beginning of file;
                %            round((20-cnt_skp)*samplingFreq*1e-3) is the number of samples to skip, in order to start the phase tracking at the beginning of the next data bit (valid only for GPS).
        end
        
        % Get a vector with the C/A code sampled 1x/chip
        %caCode = generateCAcode(channel(channelNr).PRN);
        caCode=GenerateLocCode(channel(channelNr).PRN,settings.GNSS_signal);
        % Then make it possible to do early and late versions
        caCode = [caCode(end) caCode  caCode(1)];
        
        %--- Perform various initializations ------------------------------
        
        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;
        
        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        
        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;
        % LAF loop counter
        LAFloopCnt=1;
        % CN0 loop counter
        CN0loopCnt=1;
        % Flag initialization to verify if the single channel is locked
        ChannelLocked=true;
        %=== Process the number of specified code periods =================
        loopCnt=0;
        % for loopCnt =  1:codePeriods
        while (loopCnt < codePeriods) && (ChannelLocked==true)
            loopCnt=loopCnt+1;
            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 200) == 0)
                try
                    waitbar(loopCnt/codePeriods, ...
                        hwb, ...
                        ['Tracking: Ch ', int2str(channelNr), ...
                        ' of ', int2str(settings.numberOfChannels), ...
                        '; PRN#', int2str(channel(channelNr).PRN), ...
                        '; Completed ',int2str(loopCnt), ...
                        ' of ', int2str(codePeriods), ' codePeriods']);
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end
            
            %% Read next block of data ------------------------------------------------
            % Find the size of a "block" or code period in whole samples
            
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;
            
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
            
            % Read in the appropriate number of samples to process this
            % interation
            %             [rawSignal, samplesRead] = fread(fid, ...
            %                 blksize, settings.dataType);
            %             rawSignal = rawSignal';  %transpose vector
            
            [rawSignal,samplesRead] = fread(fid,settings.NSample*blksize,settings.dataType);
            
            if strcmp(settings.signalType,'IQ')==1 %(settings.IF==0) % IQ sampling
                %I = rawSignal(1:2:settings.NSample*blksize-1)';
                %Q = rawSignal(2:2:settings.NSample*blksize)';
                rawSignal = rawSignal(1:2:settings.NSample*blksize-1) + 1i*rawSignal(2:2:settings.NSample*blksize); %%% EF:  complex representation of data
                rawSignal= rawSignal.';


                %if did not read in enough samples, then could be out of data - better exit
                if (samplesRead ~= settings.NSample*blksize)
                    disp('Tracking: Not able to read the specified number of samples, exiting...')
                    %fclose(fid);
                    return
                end
            else % IF sampling
                %[rawSignal,samplesRead] = fread(fid,blksize,settings.dataType);
                rawSignal=rawSignal';
                % If did not read in enough samples, then could be out of
                % data - better exit
                if (samplesRead ~= settings.NSample*blksize)
                    disp('Tracking: Not able to read the specified number of samples  for tracking, exiting!')
                    %fclose(fid);
                    return
                end
            end
            
            %% Set up all the code phase tracking information -------------------------
            % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);
            
            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);
            
            remcodephase_old=remCodePhase;
            remCodePhase = (tcode(blksize) + codePhaseStep) - settings.codeLength;
            
            %% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            
            %% Accumulated Phase measurements for Carrier Smoothing
            %inst_phase = 2*pi*carrFreq*Tc;
            if (loopCnt>1)
                trigarg_doppler = ((carrFreq-settings.IF)*2*pi.*time) + trackResults(channelNr).rem_fraction(loopCnt-1);
                trackResults(channelNr).rem_fraction(loopCnt) = rem(trigarg_doppler(blksize+1),(2 * pi));
                trackResults(channelNr).rem_integer(loopCnt)  = trackResults(channelNr).rem_integer(loopCnt-1) + fix(trigarg_doppler(blksize+1)/(2*pi));
                %accumulated_phase(loopCnt) =  trackResults(channelNr).rem_integer(loopCnt)+trackResults(channelNr).rem_fraction(loopCnt)/(2*pi);
                %trackResults(channelNr).accumulated_phase_scint(loopCnt) = trackResults(channelNr).accumulated_phase_scint(loopCnt-1)+inst_phase;
            else
                trigarg_doppler = (carrFreq-settings.IF)*2*pi.*time;
                trackResults(channelNr).rem_fraction(loopCnt) = rem(trigarg_doppler(blksize+1),(2 * pi));
                trackResults(channelNr).rem_integer(loopCnt)  = fix(trigarg_doppler(blksize+1)/(2*pi));
                %accumulated_phase(loopCnt) = trackResults(channelNr).rem_integer(loopCnt)+trackResults(channelNr).rem_fraction(loopCnt)/(2*pi);
                %trackResults(channelNr).accumulated_phase_scint(loopCnt) = inst_phase;
            end
            %%
            
            % Finally compute the signal to mix the collected data to bandband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));
            
            carrExp = carrSin + 1i*carrCos;

            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            BasebandSignal = carrExp .* rawSignal;  
            qBasebandSignal = imag(BasebandSignal);
            iBasebandSignal = real(BasebandSignal);
%             if settings.IF==0 % Baseband signal
%                 iBasebandSignal = I.*carrSin - Q.*carrCos;
%                 qBasebandSignal = I.*carrCos + Q.*carrSin;
%             else
%                 qBasebandSignal = carrCos .* rawSignal;
%                 iBasebandSignal = carrSin .* rawSignal;
%             end
            
            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
            %% Find PLL error and update carrier NCO ----------------------------------
            
            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);
            
            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;
            
            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;
            
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;
            
            %% Find DLL error and update code NCO -------------------------------------
            %codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) /(sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            
            codeError =0.5*((( I_E^2 + Q_E^2 ) - ( I_L^2 + Q_L^2 ))/(( I_E^2 + Q_E^2 ) + (I_L^2 + Q_L^2 )));
            
            
            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;
            
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;
            
            %% LAF part ===================================================
            if settings.enableLAF==true
                % Multicorrelator or MulticorrelatorABS or
                % parallel computing version Multicorrelator_par or
                % MulticorrelatorABS_par
                [new_u, new_d]=Multicorrelator(caCode(2:end-1),qBasebandSignal,iBasebandSignal,codePhaseStep,remcodephase_old,...
                    blksize,minlag,stepp,maxlag);
                
                %tau=[-(Nallignment-1)/2:(Nallignment-1)/2]*step/codeFreq;
                %%%%%%tau=[-(length(new_d)-1)/2:(length(new_d)-1)/2]'*step/codefreq;
                %tau_m(:,loopCnt)=tau;
                %%%
                if rem(loopCnt,settingsMPDD.mvgAvgtime)==0
                    
                    % add new correlation calculation
                    u1=(u1+new_u)/settingsMPDD.mvgAvgtime;
                    d1=(d1+new_d)/settingsMPDD.mvgAvgtime;
                    
                    
                    %% Alignment
                    if alignment==true
                        [~, indexmax]=max(u1);
                        u=MaxAlignment(u1,Nallignment,indexmax);
                        %[~, indexmax]=max(d1);
                        d=MaxAlignment(d1,Nallignment,indexmax);
                    else
                        u=u1;
                        d=d1;
                    end
                    
                    %% Normalization
                    %u=u/max(u);
                    %d=d/max(d);
                    %% LAF procedure
                    [LAF]=LAF_SQM(u,d,M,version);
                    trackResults(channelNr).LAF(LAFloopCnt) = LAF;
                    LAFloopCnt=LAFloopCnt+1;
                    
                    %% azzeramento ogni mvgAvg time?
                    u1=zeros(N,1); %u(ceil(N/2))=1;
                    d1=zeros(N,1); %d(ceil(N/2))=1;
                    
                else
                    u1=u1+new_u;
                    d1=d1+new_d;
                end
                
            end
            
            % PLL status and CN0
            [PLL_st,IP,QP,CN0aus] = control_PLL_v2_EF(I_P,Q_P,loopCnt,35,settings.samplingFreq,settings.T_int,settingsMPDD.CNo_WINDOW,...
                FIFO_IP(channelNr,:),FIFO_QP(channelNr,:));
            FIFO_IP(channelNr,:)=IP; FIFO_QP(channelNr,:)=QP;
            
            % Saving CN0 values...
            if CN0aus.CN0ready==true
                trackResults(channelNr).CN0.CNo_Emanuela(CN0loopCnt) = CN0aus.CNo_Emanuela;
                trackResults(channelNr).CN0.CNo_bluebook(CN0loopCnt) = CN0aus.CNo_bluebook;
                trackResults(channelNr).CN0.CNo_SNV(CN0loopCnt) = CN0aus.CNo_SNV;
                trackResults(channelNr).CN0.CNo_MM(CN0loopCnt) = CN0aus.CNo_MM;
                trackResults(channelNr).CN0.CNo_Bea(CN0loopCnt) = CN0aus.CNo_Bea;
                CN0loopCnt=CN0loopCnt+1;
            end
            
            if (PLL_st ~= 0)
                disp('PLL not locked, C/No below the threshold, the channel might be reassigned to another satellite, or try with a new fast acquisition...')
                ChannelLocked=false;%break;
                channel(channelNr).status='-'; 
            end;
            %==============================================================
            %% Record various measures to show in postprocessing ----------------------
            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid)/settings.BytePerSample/settings.NSample;
            
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;
            trackResults(channelNr).remCodePhase(loopCnt)   = remCodePhase;
            trackResults(channelNr).codeFreq(loopCnt)       = codeFreq;
            
            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
        end % for loopCnt
        
        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;
        
    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)

end
