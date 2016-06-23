function [transmitTime_increment] = Compute_SatTX_time_GPS(channel, settings)

%Author Gianluca Falco & Marco Rao
%Version 1.0
% Compute_SatTX_time_GPS: at a beginning of a Subframe we know that the
%                         message has been broadcast by the satellite at TOW. Since we want to
%                         compute the Satellite position at the clock time, we have to update the
%                         new time of trasmission that will be TOW+ delay (delay=time between beginning of the
%                         subframe and the clock time computed in) [s]
% Input parameters:  - GPS channels
%                    
% Output parameters: - new transmit time [s]

num_sat                 = length(channel);
transmitTime_increment  = zeros(num_sat,1);

Tc = settings.T_int;
f0 = settings.f0_L1;

for n = 1:num_sat
    doppler= channel(n).carrFreq(channel(n).Prec_clock_bit) - settings.IF;
    transmitTime_increment(n) = ...
        channel(n).eph.TOW + ...
        channel(n).bit_startpage_2clock * Tc + ...
        channel(n).AbsoluteSample_residual_samples/settings.samplingFreq + ...
        channel(n).codePhase_residual/1.023e6+...
        channel(n).AbsoluteSample_residual_samples/settings.samplingFreq*doppler/f0; % computed in sec;
end
