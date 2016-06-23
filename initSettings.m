function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure).

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 30*60*1000;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 7;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
[filename, pathname] = uigetfile({'*.bin;*.grab','"DONE" sample files (*.bin;*.grab)'}, 'Select a file',...
    'MultiSelect', 'on');
if ischar(filename)==0 % Only 1 dataset or more?
    settings.N_of_datafile      = size(filename,2);
    filename                    = fullfile(pathname, filename);
    settings.fileName           = cell(1,settings.N_of_datafile);
    settings.fileName           = {filename};
else
    settings.N_of_datafile      = 1;
    filename                    = fullfile(pathname, filename);
    settings.fileName{1,1}      = cell(1,1);
    settings.fileName{1,1}      = filename;
end

% Front-End settings
% feID =
% 1 - SiGe
% 2 - Hieu
% 3 - Thuan
% 4 - USRP
% 5 - TEXBAT
% 6 - Nfuels IQ
% 7 - Stanford USRP 1
% 8 - Stanford USRP 2
feID=8;
settings=FrontEndsettings(settings,feID);

%% Constellation characteristics
% GNSS_signal = 'GPS_L1';  % For processing GPS L1 C/A signals
% GNSS_signal = 'GAL_E1b'; % For processing Galileo E1b (OS) signals
% GNSS_signal = 'GAL_E5a';  % For processing Galileo E5a (OS) signals
% GNSS_signal = 'GAL_E5b';  % For processing Galileo E5a (OS) signals
settings.GNSS_signal = 'GPS_L1';

% Setting of Integration time [s] for the tracking (NON cambiarlo adesso!)
switch settings.GNSS_signal
    case ('GPS_L1')
        settings.T_int = 1e-3; % [s] coherent integration time
        settings.codeFreqBasis      = 1.023e6;      %[Hz]
        % Define number of chips in a code period
        settings.codeLength         = 1023;
    case ('GAL_E1b')
        settings.T_int = 4e-3; % [s] coherent integration time = 1 primary code period (= 4 ms)
        settings.codeFreqBasis      = 1.023e6;      %[Hz]
        % Define number of chips in a code period
        settings.codeLength         = 4092;
    case ('GAL_E5a')
        settings.T_int = 1e-3; % [s] coherent integration time = 1 primary code period (= 1 ms)
    case ('GAL_E5b')
        settings.T_int = 1e-3; % [s] coherent integration time = 1 primary code period (= 1 ms)
    otherwise
        error('err:unkSignal','ERROR: Unknown GNSS signal!');
end


%%
% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only.
settings.seek_sec=73*60;  %[s];
%if settings.IF==0
settings.skipNumberOfBytes     = ceil(settings.seek_sec*settings.samplingFreq*settings.NSample*settings.BytePerSample);
%else
%    settings.skipNumberOfBytes     = ceil(settings.seek_sec*settings.samplingFreq); %1000*16.368e6;%3*16.368e6;
%end



%% Acquisition settings ===================================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
settings.acqSatelliteList   = [ 8  10 11  14  22 18 32];%[1 4 11 15 18 21 24];%[2 6 12 17 19 24 25];%   [2 6 12 17 2 19];%    %[2 5 7:12 14:20 22 24:32];%[2 6 12 17 19 24 25];%[2 5 12 20 21 25 29 31];%[1, 4, 11, 15, 18, 21, 24];%[1 4 11 15 18 21 24];%[1 4 6 7 10 11 15 18 19 20 21 24 29];%[6 4 24 21 11 1  19 10 20 29 7];%[5 7 9 13 20 28 30];%[11 18 21 22 27];%[1 4 11 16 21 22 27]; % [1,4,19,21,27]; %[1:3 5:26 28:32]; % [1:15 17 18 20 21 23:26 28:32]         %[PRN numbers]
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 14; % 14           %[kHz]
settings.acqSearchStep      = 250; % 250          %[Hz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 1.2;%1.2;% 1.8; %2.0;
settings.Non_Coh_Sums       = 5;

%% FLL settings ===========================================================
settings.fll_time                = 5e3; % [ms] for the FLL
settings.fllGain                 = 10;
settings.fllNoiseBandwidth       = 25;  % 25    [Hz]

%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;       %[Hz]
settings.dllCorrelatorSpacing    = .5;%0.5;     %[chips]

% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 25; %25      [Hz]

%% Navigation solution settings ===========================================

% Period for calculating pseudoranges and position
settings.navSolPeriod       = 1000;         %[ms]

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 10;           %[degrees 0 - 90]
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off | 1 - On
% Enable/dissable use of ionospheric correction
settings.useIonoCorr        = 1;            % 0 - Off | 1 - On
settings.searchStartOffset  = 0*1e3; % Offset for function FindPreambles (usually =0). CHANGE THIS PARAMETER IF YOU DON'T FIND THE PREAMBLES!!!!

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .

% Reference position for Thuan' simulator
%Ref_User_Pos=[4472571.24804116,601624.111174408,4492080.98216451];
% ISMB Reference position
Ref_User_Pos=[4472420.8024 601433.9839 4492695.7516];  %TRUE position

Ref_User_Pos=1.0e+06*[-2.700404441100000; -4.292605170500000;3.855137392300000]; % TRUE position of the roof antenna in STANFORD

[lat,lon,~]=cart2geo(Ref_User_Pos(1),Ref_User_Pos(2),Ref_User_Pos(3),5); % llh=[45.060256 7.661107 0]
utmzone=findUtmZone(lat,lon);
[E,N,U]=cart2utm(Ref_User_Pos(1),Ref_User_Pos(2),Ref_User_Pos(3),utmzone); % Reference position
%E=NaN; N=E; U=N;
settings.truePosition.E     = E;%nan;
settings.truePosition.N     = N;%nan;
settings.truePosition.U     = U;%nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 0;            % 0 - Off
% 1 - On

%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time
settings.f0_L1              = 1575.42e6;            % [Hz] L1 radio-frequency
settings.lambdaL1           = settings.c/settings.f0_L1; % [m] wavelength of L1

%% RAIM settings ==========================================================
% RAIM scheme, probabilities (i.e. %  pfa = 6.6667e-5; pmd = 0.001;)
settings.RAIM.type=0; % 0-RS (Residual method), 1-SS (Solution separation method)
settings.RAIM.enableRAIM=true;
settings.RAIM.PfaLocal=1/30000;
settings.RAIM.Pmd=0.001;%5e-5;%1e-3*10;
settings.RAIM.HAL=50; % [m]
settings.RAIM.RecomputePVT=5;

% Computation of the number of standard deviations corresponding to the specified Pmd
k=0;
Pmd=1;
while Pmd>settings.RAIM.Pmd
    k=k+1;
    Pd=normcdf(k,0,1)-normcdf(-k,0,1);
    Pmd=1-Pd;
end
settings.RAIM.N_of_sd=k;
%%%%%%%%%%settings.RAIM.X=norminv(settings.RAIM.Pmd,0,1);
%%%%%%%%%%settings.RAIM.l1=sqrt(20.5);settings.RAIM.l2=sqrt(30.5); % 555.6;  % HAL for  NPA  (in meters)
%% LAF settings ==========================================================
settings.enableLAF=true;
%% DOP control
settings.DOPcontrol=false;
settings.GDOPbudget=.05;%9; Maximum geometric error allowable

%% Carrier Smoothing settings ==========================================================
settings.CS.enable_smoothL1    = 0;
settings.CS.smoothL1_weights   = 100;
end
