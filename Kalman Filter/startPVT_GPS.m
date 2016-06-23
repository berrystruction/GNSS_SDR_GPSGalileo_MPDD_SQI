function [startClock_PVT,channels_GPS]=startPVT_GPS(channels_GPS,num_samples_Rate,settings)
% Author Gianluca Falco & Marco Rao
% Update Nicola Linty, works both for L1 1 ms, L1 20 ms, L2 20 ms
% Version 1.0
% StartPVT_GPS     allows to compute the first clock instant to perform a PVT
%                  solution..This common clock is unique for all the channels
% Input parameters:  - GPS channels
%                    - nominal # of samples at the specific Rate
% Output parameters: - Calculated Clock time [as # of samples]
%                    - GPS channels   

NoChannels = length(channels_GPS);

for ij=1:NoChannels   
    [pos] = channels_GPS(ij).WordStart(1); % fist Subframe counted as # of bits
    if (settings.T_int == 0.001)
        [pos_bit_NextSub_frame] = channels_GPS(ij).WordStart(1)+ 300*20-1; % next Subframe start
    elseif (settings.T_int == 0.02)
        [pos_bit_NextSub_frame] = channels_GPS(ij).WordStart(1)+ 300-1; % next Subframe start
    end
    channels_GPS(ij).Word_pos = 1;
    idx = 1;
    
    while(num_samples_Rate*idx<channels_GPS(ij).absoluteSample(pos))
        idx = idx+1; % increment the clock till we are inside the first Subframe
    end
    %check if the clock is inside the first SubFrame
    if(num_samples_Rate*idx<=channels_GPS(ij).absoluteSample(pos_bit_NextSub_frame))
        startClock_PVT = num_samples_Rate*idx;
    else
        disp('Error in the Page counter..Exiting')
    end
end
             
return;