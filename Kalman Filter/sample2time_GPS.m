function [timeDifference,travelTime,pseudorange,delta_pseudorange] = sample2time_GPS(channel,Doppler,startOffset, settings)

%Author Gianluca Falco & Marco Rao
%Version 1.0
% sample2time_GPS: allows to compute the pseudoranges [m],delta_pseudoranges [m], traveltime [ms] 
%                  and the timeDifference between the different Satellites' channels [s]
%                      
% Input parameters:  - GPS channels
%                    - Doppler freq [Hz]
%                    - startOffset is the min_travelTime in ms
% Output parameters: - pseudoranges [m]
%                    - delta_pseudoranges [m]
%                    - traveltime [ms]
%                    - timeDifference between between the different Satellites' channels [s]

f_sampling = settings.samplingFreq;
Tc = settings.T_int;
lambda = settings.lambdaL1;
f0 = settings.f0_L1;

num_satellite   = length(channel);
timeDifference  = zeros(num_satellite,1);
travelTime      = zeros(num_satellite,1);
pseudorange     = zeros(num_satellite,1);
delta_pseudorange = zeros(num_satellite,1);

for ind = 1:numel(channel)
    time_distance(ind) = channel(ind).AbsoluteSample_residual_samples / f_sampling + ...
        channel(ind).bit_startpage_2clock * Tc+...
        channel(ind).codePhase_residual/settings.codeFreqBasis + ...
        channel(ind).AbsoluteSample_residual_samples*Doppler(ind) / f_sampling / f0; % This last line propagates the Doppler
end

[~, pos] = max(time_distance);
%clear dump;
% Time difference with respect to the maximum travel time
for ind = 1:numel(channel)
    timeDifference(ind) = (channel(ind).AbsoluteSample_residual_samples - channel(pos).AbsoluteSample_residual_samples) / f_sampling + ...
        (channel(ind).bit_startpage_2clock - channel(pos).bit_startpage_2clock) * Tc + ...
        (channel(ind).codePhase_residual - channel(pos).codePhase_residual)/settings.codeFreqBasis + ...
        (channel(ind).AbsoluteSample_residual_samples*Doppler(ind) - channel(pos).AbsoluteSample_residual_samples * Doppler(pos)) / f_sampling / f0; % This last line propagates the Doppler
end
    
[~, ind_closer] = min(abs(timeDifference)); %find the minum travel_time = the closest satellite  to the user
%clear dump;
increments = timeDifference - timeDifference(ind_closer);

for ind = 1:numel(channel)
%% Thuan - modified
%     if (settings.signal == 'GPS L2')
%         if (abs(increments(ind))*1e3 > 20) % This can happen if the first bit of the navigation message is discarded to du Viterbi fails in decoding
%             increments(ind) = rem(increments(ind)*1e3,20)/1e3;
%         end
%     end
%% end modified
    % Code_Index_Start = channel(ind).Prec_clock_bit;
    travelTime(ind) = startOffset + abs(increments(ind)) * 1000;%  ...
%%
    %  - channel(ind).codePhase(Code_Index_Start)/1.023e6*1000; % old
    %  version, corrected by Thuan
    pseudorange(ind) = travelTime(ind) * settings.c / 1000; 
    delta_pseudorange(ind) = Doppler(ind) * lambda;
end
% 'c'