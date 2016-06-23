function [channels_GPS,Doppler] = Counter_clock_Update_GPS(channels_GPS, clock, settings, sol_index, Drift)

%Author Gianluca Falco & Marco Rao
%Version 1.0
% Counter_clock_Update_GPS: 1)allows to compute the time difference (measured
%                             in samples) between the clock time and the
%                             beginning of the GPS subframe
%                           2) Doppler computation
% Input parameters:  - GPS channels
%                    - clock [as # of samples]
% Output parameters: - Calculated Doppler freq [Hz]
%                    - GPS channels where the following variables have been added to the previous channels struct:
%                         .AbsoluteSample_residual_samples: measurements of the time distance between the clock instant and the nearest
%                                                           data_bit transition [# samples]
%
%                         .bit_startpage_2clock: difference in bit between the beginning of a GPS subframe
%                                                 and the nearest data bit to the clock time [# bits]
%
%                         .Samples_Startword_2clock: number of samples between the clock time and the beginning
%                                                    of a GPS subframe [#
%                                                    samples]

num_sat = length(channels_GPS);

for ij=1:num_sat
    %detection of the residual samples between the beginning of the nearest previous data bit
    %and the clock time sample
    while channels_GPS(ij).ind<length(channels_GPS(1).absoluteSample) && channels_GPS(ij).absoluteSample(channels_GPS(ij).ind)-clock < 0
        channels_GPS(ij).ind = channels_GPS(ij).ind+1;
    end
    [pos_Prec_bit] = channels_GPS(ij).ind-1;
    [pos_Next_bit] = channels_GPS(ij).ind;
    clear i_tmp;
    
    %arrange the two residual samples
    if (channels_GPS(ij).absoluteSample(pos_Prec_bit) > channels_GPS(ij).absoluteSample(pos_Next_bit))
        a            = pos_Prec_bit;
        pos_Prec_bit = pos_Next_bit;
        pos_Next_bit = a;
    end
    channels_GPS(ij).Prec_clock_bit = pos_Prec_bit;
    channels_GPS(ij).Next_clock_bit = pos_Next_bit;
    % compute the residual samples
    channels_GPS(ij).AbsoluteSample_residual_samples = clock - channels_GPS(ij).absoluteSample(pos_Prec_bit) + 1; % measured in samples
    
    % compute the related code-offset
    channels_GPS(ij).codePhase_residual = channels_GPS(ij).remCodePhase(pos_Prec_bit);
    
    % Now we look for the bit position
    bit_SubFrame_start = channels_GPS(ij).WordStart(channels_GPS(ij).Word_pos); % measured in bits
    % compute the distance in bit between the beginning of the subframe and
    % the closest data_bit(pos_Prec_bit) to the clock time
    bit_expiredfrom_SubFrame = pos_Prec_bit - bit_SubFrame_start + 1; %measured in bits
    
    channels_GPS(ij).bit_startpage_2clock = bit_expiredfrom_SubFrame; % measurement of time distance between the clock instant and the beginning of the Subframe in bits
    
    % compute the samples (time) difference between the nearest bit and the
    % beginning of the subframe
    diff_bit_wordstart = channels_GPS(ij).absoluteSample(pos_Prec_bit) - ...
        channels_GPS(ij).absoluteSample(channels_GPS(ij).WordStart(channels_GPS(ij).Word_pos));
    
    %adding the residual bits
    channels_GPS(ij).Samples_Startword_2clock = diff_bit_wordstart+channels_GPS(ij).AbsoluteSample_residual_samples;
    
    % Doppler measurements
    Doppler(ij) = channels_GPS(ij).carrFreq(pos_Prec_bit) - settings.IF;
    
    % Drift correction
    Doppler(ij) = Doppler(ij)*(1-(Drift/settings.c));
    
    if (settings.CS.enable_smoothL1)
        time = (0:channels_GPS(ij).AbsoluteSample_residual_samples) ./ settings.samplingFreq;
        trigarg_doppler = (Doppler(ij) * 2.0 * pi) .* time + channels_GPS(ij).rem_fraction(pos_Prec_bit);
        channels_GPS(ij).phase_TR_fract = rem(trigarg_doppler(end), (2*pi))*settings.lambdaL1/(2*pi);
        channels_GPS(ij).phase_TR_fract_rad = rem(trigarg_doppler(end), (2*pi));
        channels_GPS(ij).phase_TR_int = channels_GPS(ij).rem_integer(pos_Prec_bit)+fix(trigarg_doppler(end)/(2*pi));
        channels_GPS(ij).Lpc(1,sol_index) = channels_GPS(ij).phase_TR_fract + settings.lambdaL1*channels_GPS(ij).phase_TR_int;     %meters
    end
end
