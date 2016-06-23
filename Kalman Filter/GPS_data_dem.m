function [channels, subFrameStart] = GPS_data_dem(channels,Tc,settings)
% CHECK IOD CAREFULLY
NoChannels = length(channels);
[subFrameStart] = GPS_findPreambles(channels,Tc,settings);

%% Decode ephemerides =====================================================

for channelNr = 1:NoChannels
    channels(channelNr).WordStart = subFrameStart(channelNr);
    
    %=== Convert tracking output to navigation bits =======================
    if (Tc == 0.001)
        %--- Copy 5 sub-frames long record from tracking output ---------------
        navBitsSamples = channels(channelNr).I_P(subFrameStart(channelNr) - 20 : ...
            subFrameStart(channelNr) + (300 * 5 * 20) -1)'; %%where 20 is the number of samples in 20ms
        
        
        %--- Group every 20 vales of bits into columns ------------------------
        navBitsSamples = reshape(navBitsSamples, ...
            20, (size(navBitsSamples, 1) / 20));
        
        %--- Sum all samples in the bits to get the best estimate -------------
        navBits = sum(navBitsSamples);
    elseif (Tc == 0.02)
        %--- Copy 5 sub-frames long record from tracking output ---------------
        navBits = channels(channelNr).I_P(subFrameStart(channelNr) - 1 : ...
            subFrameStart(channelNr) + (300 * 5 ) -1)'; %%where 1 is the number of samples in 20ms
    end
    
    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0);

    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    
    %=== Decode ephemerides and TOW of the first sub-frame ================
    [channels(channelNr).eph, TOW] = ...
                            ephemeris(navBitsBin(2:1501)', navBitsBin(1));
    channels(channelNr).eph.TOW=TOW; % MODIFICA CHE HO FATTO IO MATTIA!!  
    
% ephemeris2 contains demodulation for ionospheric parameters                        
%       [channels(channelNr).eph, TOW] = ...
%                             ephemeris2(navBitsBin(2:1500*25+1)', navBitsBin(1));                    

end


return