function [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_Gal_E1c_BOC(SatVNumber,MemoryMsg)
% GenModulation_Gal_E1c_BOC
%   
%   One period of the Galileo E1 signal (100 ms) of the specified satellite is
%   generated around baseband(BB) and with the minimum required sampling
%   frequency (2.046 MHz). 
%   Specifications:
%   - BOC(1,1) modulation;
%   - primary and secondary code on C channel;
% -------------------------------------------------------------------------
% Inputs:
% SatVNumber        Desired satellite
% -------------------------------------------------------------------------
% Outputs:
%   CodeBaseCell        cell array containing:
%       NumChannels     number of code channels in the modulation
%       {CodeBasePeriod}  basic period of each code channel (Vector containes the samples of 1 period of the GPS L1
%                   C/A code)
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
% B. Motella
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

FreqCode = 2*1.023e6;                               % minimum required sampling frequency

load('Gal_E1_Codes.dat', '-mat');                   % matrix contains the satellites PRN codes loading 

Cp = pr_E1_C(SatVNumber,:);                         % primary code - channel C
Cs = sc_E1_C;                                       % secondary code - channel C
Cp = repmat(Cp,1,length(Cs));                       % primary code - repetition to fit the secondary code
Cs = sc_E1_C(floor((0:length(Cs)*4092-1)/4092)+1);  % secondary code - resampling
Cps = Cp.*Cs;                                       % code channel C

% Power normalization
Pcode = mean(abs(Cps).^2); % Assuming a theoretical mean value equal to 0
Cps = Cps/sqrt(Pcode);

% introduction of the BOC(1,1)
CodeBasePeriod = zeros(1,4092*2*25);
CodeBasePeriod(1:2:end) = -Cps;               % pilot channel only 
CodeBasePeriod(2:2:end) = Cps;              % pilot channel only 

if MemoryMsg.DataFlag == 1
    MemoryMsg.BitDuration = 1/0;
    MemoryMsg.AssignedBitCounter = 1;
    MemoryMsg.CurrentBitSign = MemoryMsg.NavMessage(1);
end

CodeBaseCell = {1  CodeBasePeriod  100e-3};
