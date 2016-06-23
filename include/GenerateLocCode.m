function [Code, Rc] = GenerateLocCode(PRN,GNSS_signal)
%
% GenerateLocCode: this function generates the local code for the
% acquisition, depending on the chosen GNSS signal (band & channel) and PRN. 
% It can be customized in order to call user refined code genrators and
% return the code rate (Rc).

%global GNSS_signal;

switch GNSS_signal
    case ('GPS_L1')        
        % GPS CA code generator
        Code=CAgen(PRN);        
%         warning('4 GPS C/A code periods (4 ms)!')
%         Code = [Code Code Code Code]; % 4 ms
        Rc = 1.023*1e6;              % Nominal GPS C/A code rate.
        
    case ('GAL_E1b')
        
        if (PRN == 51)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Galileo GIOVE-A E1b code generator
            MemoryMsg.DataFlag = 0; % No data
            [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_GioveA_E1b_BOC(PRN,MemoryMsg);
            Code=CodeBaseCell{2};       % One code period (4 ms)
            %     [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_GioveA_E1c_BOC(PRN,MemoryMsg);
            %     Code=CodeBaseCell{2};       % One code period (200 ms)
            %     Code = Code(1:8184*2);      % One primary code period (8 ms)
            Rc = FreqCode;              % Nominal code rate.

        elseif (PRN == 52)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Galileo GIOVE-B E1b code generator
            MemoryMsg.DataFlag = 0; % No data
            [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_GioveB_E1b_BOC(PRN,MemoryMsg);
            Code=CodeBaseCell{2};       % One code period (4 ms)
            %     [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_GioveB_E1c_BOC(PRN,MemoryMsg);
            %     Code=CodeBaseCell{2};       % One code period (200 ms)
            %     Code = Code(1:8184*2);      % One primary code period (8 ms)
            Rc = FreqCode;              % Nominal code rate.

        else
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Galileo E1b code generator
            MemoryMsg.DataFlag = 0; % No data
            [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_Gal_E1b_BOC(PRN,MemoryMsg);
            Code=CodeBaseCell{2};       % One code period (4 ms)
            %     [CodeBaseCell, FreqCode, MemoryMsg] = GenModulation_Gal_E1c_BOC(PRN,MemoryMsg);
            %     Code=CodeBaseCell{2};       % One code period (100 ms)
            %     Code = Code(1:8184);        % One primary code period (4 ms)
            Rc = FreqCode;              % Nominal code rate.
            
        end

        
    case ('GAL_E5a')           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Galileo E5a-I code generator   
        Rc = 10.23e6; % Nominal code rate.
        T_code = 1e-3; % [s] code length
        N = Rc*T_code;  % 10230 % Number of generated code chips.
        % Data channel (no data modulation is considered, 1 primary code period)  
        [Code] = E5_code_gen(PRN, 'E5aI', N); % Galileo E5 code generator            

        
    case ('GAL_E5b')           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Galileo E5b-I code generator   
        Rc = 10.23e6; % Nominal code rate.
        T_code = 1e-3; % [s] code length
        N = Rc*T_code;  % 10230 % Number of generated code chips.
        % Data channel (no data modulation is considered, 1 primary code period)  
        [Code] = E5_code_gen(PRN, 'E5bI', N); % Galileo E5 code generator       
        
        
    otherwise
        error('err:unkSignal','ERROR: Unknown GNSS signal!')
end





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GIOVE-A/B E5 code generator
% Rc = 10.23e6; % Mcps
% T_code = 1e-3; % [s] code length
% N = Rc*T_code;  % 10230 % Number of generated code chips.
% 
% [Code] = E5_code_gen_GIOVE(PRN, 'E5bQ', N); % GIOVE-A/B E5 code generator



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GPS L2C CM (20 ms) code generation v1
% Rc = 1.023*1e6;              % Nominal L2C combined CM/CL code rate.
% 
% L_CodeCM = 10230;                                       % CM code length, 10230 chip = 20 ms
% L_CodeCL = 767250;                                      % CL code length, 767250 chip = 1.5 s 
% 
% [CodeCM] = CMCL_Generator(PRN,'CM',L_CodeCM);           % CM code generation 
% CodeCM_rep = CodeCM;
% % CodeCM_rep = repmat(CodeCM,1,L_CodeCL/L_CodeCM);      % CM code repetition to fit the CL code length
% % [CodeCL] = CMCL_Generator(PRN,'CL',L_CodeCL);         % CL code generation 
% 
% Code = zeros(1,2*length(CodeCM_rep));                   % L2C signal vector
% 
% Code(1:2:end) = CodeCM_rep;                             % CM & CL multiplexing 
% % Code(2:2:end) = CodeCL;                               % CM & CL multiplexing 
