function settings=FrontEndsettings(settings,feID)
% INPUT:
% settingsID: front end ID, integer for identifying front-end
% OUTPUT:
% settings: front-end data inside settings structure

switch feID
    case {1} % SiGe
        settings.signalType= 'IF';
        settings.NSample = 1;
        % Intermediate frequency
        settings.IF                 = 4.1304e6;  % [Hz]
        settings.dataType           = 'int8';
        % Sampling and code frequencies
        settings.samplingFreq       = 16.3676e6; % [Hz]
        
    case {2} % Hieu
        settings.signalType= 'IF';
        settings.NSample = 1;
        % Intermediate frequency
        settings.IF                 = 4.123e6; % [Hz]
        settings.dataType           = 'int8';
        % Sampling and code frequencies
        settings.samplingFreq       = 16.367e6; % [Hz]
        
    case {3} % Thuan
        settings.signalType= 'IF';
        settings.NSample = 1;
        % Intermediate frequency
        settings.IF                 = 4.1304e6;  % [Hz]
        settings.dataType           = 'int8';
        % Sampling and code frequencies
        settings.samplingFreq       = 16.368e6;  % [Hz]
        
    case {4} % USRP
        settings.signalType= 'IQ';
        settings.NSample = 2;
        % Intermediate frequency (Baseband settings)
        settings.IF                 = 0;      % [Hz]
        settings.dataType           = 'int16';
        % Sampling and code frequencies
        settings.samplingFreq       = 5e6;    % [Hz]
        
    case {5} % TEXBAT
        settings.signalType= 'IQ';
        settings.NSample = 2;
        % Intermediate frequency (Baseband settings)
        settings.IF                 = 0;      %[Hz]
        settings.dataType           = 'int16';
        % Sampling and code frequencies
        settings.samplingFreq       = 25e6;      % [Hz]
        
    case {6} % Nfuels
        settings.signalType= 'IQ';
        settings.NSample = 1;
        % Intermediate frequency (Baseband settings)
        settings.IF                 = 0;      %[Hz]
        settings.dataType           = 'int8';
        % Sampling and code frequencies
        settings.samplingFreq       = 8367600;      % [Hz]
        
    case {7} % Stanford USRP 1
        settings.signalType= 'IQ';
        settings.NSample = 2;
        % Intermediate frequency (Baseband settings)
        settings.IF                 = 0.42e6;      %[Hz]
        settings.dataType           = 'short';
        %settings.BytePerSample = 2;
        % Sampling and code frequencies
        settings.samplingFreq       = 10e6;      % [Hz]
        
    case {8} % Stanford USRP 2
        settings.signalType= 'IQ';
        settings.NSample = 2;
        % Intermediate frequency (Baseband settings)
        settings.IF                 = 0;      %[Hz]
        settings.dataType           = 'short';
        %settings.BytePerSample = 2;
        % Sampling and code frequencies
        settings.samplingFreq       = 8e6;      % [Hz]
end


% Data type used to store one sample
switch settings.dataType
    case {'int8','char','schar'}
        settings.BytePerSample = 1;
    case {'int16','short'}
        settings.BytePerSample = 2;
end


end

