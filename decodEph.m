function [eph, TOW, newactiveChnList]=decodEph(trackResults,activeChnList,subFrameStart,IONOon)
% Decode ephemerides annd ionospheric parameters
first_5_subframes=zeros(length(activeChnList),length([subFrameStart(1) - 20 :subFrameStart(1) + (1500 * 20)-1]));

for channelNr = activeChnList
    
    %=== Convert tracking output to navigation bits =======================
    % first 5 sub-frames
    first_5_subframes(channelNr,:)=[subFrameStart(channelNr) - 20 :subFrameStart(channelNr) + (1500 * 20)-1];
    %--- Copy 5 sub-frames long record from tracking output ---------------
    navBitsSamples = trackResults(channelNr).I_P(first_5_subframes(channelNr,:))';
    
    %--- Group every 20 values of bits into columns ------------------------
    navBitsSamples = reshape(navBitsSamples, ...
        20, (size(navBitsSamples, 1) / 20));
    
    %--- Sum all samples in the bits to get the best estimate -------------
    navBits = sum(navBitsSamples);
    
    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0);
    
    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    
    %=== Decode ephemerides and TOW of the first sub-frame ================
    [eph(trackResults(channelNr).PRN), TOW] = ...
        ephemeris(navBitsBin(2:1501)', navBitsBin(1));
    
    
    %--- Exclude satellite if it does not have the necessary nav data -----
    if (isempty(eph(trackResults(channelNr).PRN).IODC) || ...
            isempty(eph(trackResults(channelNr).PRN).IODE_sf2) || ...
            isempty(eph(trackResults(channelNr).PRN).IODE_sf3))
        
        %--- Exclude channel from the list (from further processing) ------
        newactiveChnList = setdiff(activeChnList, channelNr);
    else
        newactiveChnList=activeChnList;
    end
end

if IONOon==true
    %% Iono corrections
    for channelNr = activeChnList
        % For each channel, I search ionospheric params
        SV_ID=0; % initialization
        kk=0;
        while SV_ID~=56 % if we are at page 18 compute the ionosphere parameters
            kk=kk+1; % index multiplier
            
            %=== Convert tracking output to navigation bits =======================
            %--- Copy 5 sub-frames long record from tracking output ---------------
            navBitsSamples = trackResults(channelNr).I_P(first_5_subframes(channelNr,:)+kk*(1500 * 20))';%first_5_subframes(channelNr,1500))';
            
            %--- Group every 20 vales of bits into columns ------------------------
            navBitsSamples = reshape(navBitsSamples,20, (size(navBitsSamples, 1) / 20));
            
            navBits = sum(navBitsSamples);
            navBits = (navBits > 0);
            navBitsBin = dec2bin(navBits);
            
            %=== Decode ephemerides and TOW of the first sub-frame ================
            [ephIono, TOW_iono, SV_ID] = ...
                ephemeris(navBitsBin(2:1501)', navBitsBin(1));
        end
        % Copying Iono parameters...
        eph(trackResults(channelNr).PRN).alpha_0 = ephIono.alpha_0;
        eph(trackResults(channelNr).PRN).alpha_1 = ephIono.alpha_1;
        eph(trackResults(channelNr).PRN).alpha_2 = ephIono.alpha_2;
        eph(trackResults(channelNr).PRN).alpha_3 = ephIono.alpha_3;
        
        eph(trackResults(channelNr).PRN).beta_0 = ephIono.beta_0;
        eph(trackResults(channelNr).PRN).beta_1 = ephIono.beta_1;
        eph(trackResults(channelNr).PRN).beta_2 = ephIono.beta_2;
        eph(trackResults(channelNr).PRN).beta_3 = ephIono.beta_3;
        eph(trackResults(channelNr).PRN).iflag  = 1;
        
    end
end

end