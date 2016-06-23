function [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_Gal_E1bc_CBOC(SatVNumber, MemoryMsg);
% GenModulation_Gal_E1bc_CBOC
%  
%   One period of the Galileo CBOC signal foreseen for the E1 band is
%   generated for the specified satellite around baseband(BB) and with the
%   minimum required sampling frequency (1.023 MHz).
% -------------------------------------------------------------------------
% Inputs:
% SatVNumber        Desired satellite
% -------------------------------------------------------------------------
% Outputs:
% CodeBasePeriod    Vector containes the samples of 1 period of the Galileo
%                   E1 signal (100 ms) with the CBOC(6,1,1/11) modulation
% FreqCode          CodeBasePeriod sampling frequency
% -------------------------------------------------------------------------
% Outputs:
%   CodeBaseCell        cell array containing:
%       NumChannels     number of code channels in the modulation
%       {CodeBasePeriod}  Vector containes the samples of 1 period of the Galileo
%                   E1 signal (100 ms) with the CBOC(6,1,1/11) modulation
%       {PeriodLength}  length of one code period in seconds
%
%   FreqCode          CodeBasePeriod sampling frequency
%
%   MemoryMsg           structure contains the following fields (to manage the Navigation Message insertion)
%       DataFlag            Flag indicating the presence (1)  or absence (0) of the Navigation Message
%
%       NavMessage          [ModulationChannelNumber x n] Bit sequence of
%       the navigation message(s)
%
%       BitDuration         [ModulationChannelNumber x 1] Duration of 1 data bit for each code channel of the selected modulation [s] 
%
%       AssignedBitCounter  [ModulationChannelNumber x 1] Number of data bits already assigned in the transmitted SIS [1]
%
%       CurrentBitSign      [ModulationChannelNumber x 1] Sign of the last assigned bit [+/-1]
%       The last three fields are ABSENT if DataFlag == 0
%
% -------------------------------------------------------------------------
%
% This file is a part of the N-FUELS Signal Generation and Analysis Tool.
%
% All rights reserved, 2009.
% -------------------------------------------------------------------------
% Created by:
% D. Margaria
% 
% Last update:
% D. Margaria
% 
% Navigation Lab
% Politecnico di Torino / Istituto Superiore Mario Boella
%
%
% Version 2.2
% 14/10/2009
% -------------------------------------------------------------------------


% Flag_channel = 'B';     % Flag enabling the generation of only the data channel
% Flag_channel = 'C';     % Flag enabling the generation of only the pilot channel
Flag_channel = 'BC';      % Flag to enable the generation of both data and pilot channels


 N = 12;
 FreqCode = N*1.023e6;                              % minimum required sampling frequency

load('Gal_E1_Codes.dat', '-mat');                   % matrix contains the satellites PRN codes loading 

B = pr_E1_B(SatVNumber,:);                          % primary code - channel B
Cp = pr_E1_C(SatVNumber,:);                         % primary code - channel C
Cs = sc_E1_C;                                       % secondary code - channel C
B = repmat(B,1,length(Cs));                         % primary code - repetition to fit the secondary code
Cp = repmat(Cp,1,length(Cs));                       % primary code - repetition to fit the secondary code
Cs = sc_E1_C(floor((0:length(Cs)*4092-1)/4092)+1);  % secondary code - resampling
Cps = Cp.*Cs;                                       % code channel C

E_B = B(floor((0:length(B)*N-1)/N)+1);              % code channel B (upsampled)
E_C = Cps(floor((0:length(Cps)*N-1)/N)+1);          % code channel C (upsampled)

Seq_BOC61 = repmat([1 -1],1,6);                     % BOC(6,1) slot generation 
Seq_BOC11 = [ones(1,6) -ones(1,6)];                 % BOC(1,1) slot generation 

alpha = sqrt(10/11);
beta = sqrt(1/11);
Chip_data = [alpha*Seq_BOC11+beta*Seq_BOC61];                   % CBOC data chip generation 
Seq_data = repmat(Chip_data,1,length(E_B)/length(Chip_data));   % CBOC data sequence generation 
Chip_pilot = [alpha*Seq_BOC11-beta*Seq_BOC61];                  % CBOC pilot chip generation 
Seq_pilot = repmat(Chip_pilot,1,length(E_C)/length(Chip_pilot));% CBOC pilot sequence generation 

B_ch = (Seq_data.*E_B);   % data channel (B) only
C_ch = (Seq_pilot.*E_C);  % pilot channel (C) only

switch Flag_channel % Output signal selection
    
    case 'B'
        disp('Warning: only the data channel of Galileo E1 CBOC modulation is generated.')
        CodeBasePeriod = B_ch;   % data channel (B) only

    case 'C'
        disp('Warning: only the pilot channel of Galileo E1 CBOC modulation is generated.')
        CodeBasePeriod = C_ch;    % pilot channel (C) only
        
    otherwise
        CodeBasePeriod = (B_ch-C_ch); % data and pilot channels

end
 
if MemoryMsg.DataFlag == 0
    
    % Power normalization
    PcodeBC = mean(abs(CodeBasePeriod).^2); % Assuming a theoretical mean value equal to 0
    CodeBasePeriod = CodeBasePeriod./sqrt(PcodeBC);
    
    CodeBaseCell = {1  CodeBasePeriod  100e-3};
    
else   
    MemoryMsg.BitDuration = [1/250   1/0];
    MemoryMsg.AssignedBitCounter = [1  1];
    MemoryMsg.CurrentBitSign = [MemoryMsg.NavMessage(1, :) 1];
    %MemoryMsg.CurrentBitSign = MemoryMsg.NavMessage(1);
        
    % Power normalization
    PcodeB = mean(abs(B_ch).^2); % Assuming a theoretical mean value equal to 0
    PcodeC = mean(abs(C_ch).^2); % Assuming a theoretical mean value equal to 0
    B_ch = B_ch./sqrt(PcodeB);
    C_ch = C_ch./sqrt(PcodeC);

    CodeBaseCell = {2  (B_ch) 4e-3  (C_ch) 100e-3};

end
