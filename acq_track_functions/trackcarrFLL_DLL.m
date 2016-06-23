function [FLL_st,counter,carrfreq,codephase_out] = trackcarrFLL_DLL(file,settings,mycodephase,mycarrfreq,prn)%(file,f_sampling,prn, mycodephase, mycarrfreq, msec)%, plotme)

f_sampling=settings.samplingFreq;
msec=settings.fll_time;
T_int=settings.T_int;
%prn=settings.PRN;
plotme=0; % To enable plot of FLL results


% Code/Carrier Tracking
% M. Pini on the basis of D. Akos - UCBoulder - ASEN5519 - Spring 2005; (c) D. Akos
%
%Inputs
%  prn - the prn number of the SV to track
%  mycodephase - the starting sample number in the data file to start tracking at (from acq)
%  mycarrfreq - the starting carrier frequency for the SV to track
%  msec - the number of msec to run for (be sure there is enough data for this)
%  plotme - does a plot of the output if "1", otherwise no plot is drawn
%
%Outputs
%  FLL_st - FLL status, 0 not locked, 1 locked, -1 timeout
%  frequency_out: frequency estimated by the FLL
%  codephase_out: codephase estimated by the FLL
%  counter: is the number of milliseconds from the last data bit transaction - used to start the PLL with a PDI of 20 ms
%Requires
% m-file: cacode.m (to generate the C/A code, 1 sample/chip)
% input data file named "gpsdata.bin" of signed char type data format  "gpsdata_long.bin" of signed char type data format
%  "numsamp.bin" of short type data format

% fseek(file,0,-1);
% [gpsdata, samplesRead] = fread(file,50000,'schar');  % read approx. 3 mseconds of data
% gpsdata = gpsdata';

%skip through that data file to start at the appropriate sample (corresponding to code phase)
% fseek(file,(mycodephase-1)*16,'bof'); % DOUBLE!!!!
fseek(file,settings.NSample*settings.BytePerSample*(mycodephase-1),'bof'); % SCHAR!!

% %get a vector with the C/A code sampled 1x/chip
% ca=cacode(prn);
% %then make 1st value be the last chip and the last value the very first chip
% ca=[ca(1023) ca ca(1)];

% Local code generation
[ca, Rc] = GenerateLocCode(prn,settings.GNSS_signal);

%define number of chips in a code period
numchips=length(ca);
% numchips=length(ca)*20;
%then make 1st value be the last chip and the last value the very first chip
ca=[ca(end) ca ca(1)];
% ca=[ca(end) ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca(1)];

% -------------------------------------------------------------------------
% perform various initializations
%loop counter for times through the loop
loopcnt=1;
%define sampling frequency
fs=f_sampling;

% define initial code frequency and code frequency basis of NCO in chips/sec
% codefreq = 1.0230e6;
codefreq = Rc;
% NCO values initialization
% codefreq_basis = 1.023e6;
codefreq_basis = Rc;
%define number of chips in a code period
% numchips=1023;
% numchips=length(ca);
%define code phase (in chips) which is "left over"
remcodephase = 0.0;
%define early-late offset (in chips)
earlylate = 0.5;

%define carrier frequency which is used over whole tracking period (ok if apprx as noncoherent DLL)
carrfreq = mycarrfreq;
%carrierfreq_basis = mycarrfreq;
%define how much carrier phase is "left over" at the end of each code period to reset trig argument
remcarrphase = 0.0;
%--------------------------------------------------------------------------

% CARRIER TRACKING LOOP
zeta_carr = 1/sqrt(2);      % Butterworth
Bl_carr = settings.fllNoiseBandwidth; %15; %25;
wn_carr = Bl_carr/0.53;            % Natural Frequency;
k_carr = settings.fllGain; %10; %1;                 % Gain of the overall loop
% loop filter coefficients
t1_carr = k_carr/(wn_carr*wn_carr);
t2_carr = 2*zeta_carr/wn_carr;
% filter values initialization
t2_div_t1_carr = t2_carr/t1_carr;
olderrort_carr = 0;
oldcarrier_nco = 0;
% delta_carr = 1e-3;               % integration time
delta_carr = settings.T_int;               % integration time
% delta_carr = 20e-3;               % integration time
delta_div_t1_carr = delta_carr/t1_carr;
errort_carr = 0;
counter = 0;                     % counter is a variable used to count the number of milliseconds after a detection of the data bit transition. This variable is useful
% to pass the PLL the right number of samples to skip, in order start with a PDI of 20 ms
vec_freq_carr = [];
vec_errort_carr = [];
vec_nco_carr =[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CODE TRACKING LOOP
zeta = 1/sqrt(2);      % Butterworth
Bl=5;
wn=Bl/0.53;              % Natural Frequency;
k = 1;                 % Gain of the overall loop

% loop filter coefficients
t1 = k/(wn*wn);
t2 = 2*zeta/wn;
% filter values initialization
t2_div_t1 = t2/t1;
olderrort_code = 0;
oldcode_nco = 0;
% delta = 1e-3;               % integration time
delta = T_int;               % integration time
% delta = 20e-3;               % integration time
delta_div_t1 = delta/t1;
errort_code = 0;
codephase_out = mycodephase;

Freq_sum=0;
Old_Freq_sum=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%start a timer
%tic;
fprintf(['FLL processing for PRN ' num2str(prn) '... '])
% while less than the number of specified code periods, continue processing
% while (loopcnt <=  msec)
while (loopcnt*(T_int*1e3) <=  msec)
    % while (loopcnt*20 <=  msec)
    
    %since it can be time consuming, do a periodic display of current stage
    %also a good point to display debugging information
    %     if (rem((loopcnt),200)==0)
    if (rem((loopcnt*(T_int*1e3)),200)==0)
        %disp(['   Currently on interation:',int2str(loopcnt*(T_int*1e3)),' of:',int2str(msec)])
        fprintf('%dms ', loopcnt*(T_int*1e3));
        %     if (rem((loopcnt*20),200)==0)
        %         disp(['   Currently on interation:',int2str(loopcnt*20),' of:',int2str(msec)])
    end
    
    %update the phasestep based on code freq (variable) and sampling frequency (fixed)
    codephasestep=codefreq/fs;
    vet_freq(loopcnt) = codefreq;
    
    %find the size of a "block" or code period in whole samples
    %note - not starting from zero, but rather where you left off from last code period
    %this could be done using the CA code generator with sampling freq, code freq, and offset inputs
    %but this will be faster
    blksize=ceil((numchips-remcodephase)/codephasestep);
    
    %read in the appropriate number of samples to process this interation
    %     [rawdata,samplesRead] = fread(file,blksize,'int8');
    %     %     [rawdata,samplesRead] = fread(file,blksize,'double');
    %     rawdata=rawdata';
    %     codephase_out = codephase_out + blksize; % Update the number of processed samples. This parameter will be passed along to the PLL...
    
    [rawdata,samplesRead] = fread(file,settings.NSample*blksize,settings.dataType);
    
    if strcmp(settings.signalType,'IQ')==1 %(settings.IF==0) % IQ sampling
        %I = rawdata(1:2:2*blksize-1)';
        %Q = rawdata(2:2:2*blksize)';
        rawdata = rawdata(1:2:settings.NSample*blksize-1) + sqrt(-1)*rawdata(2:2:settings.NSample*blksize); %%% EF:  complex representation of data
        rawdata= rawdata.';
        
        %if did not read in enough samples, then could be out of data - better exit
        if (samplesRead ~= settings.NSample*blksize)
            disp('FLL: Not able to read the specified number of samples, exiting...')
            fclose(file);
            return
        end
        %         codephase_out = codephase_out + settings.NSample*blksize; % Update the number of processed samples. This parameter will be passed along to the PLL...
    else % IF sampling
        
        rawdata=rawdata';
        % If did not read in enough samples, then could be out of
        % data - better exit
        if (samplesRead ~= settings.NSample*blksize)
            disp('FLL: Not able to read the specified number of samples  for tracking, exiting!')
            fclose(file);
            return
        end
        %         codephase_out = codephase_out + blksize; % Update the number of processed samples. This parameter will be passed along to the PLL...
    end
    
    
    codephase_out = codephase_out + blksize; % Update the number of processed samples. This parameter will be passed along to the PLL...
    
    % % % %     %if did not read in enough samples, then could be out of data - better exit
    % % % %     if (samplesRead ~= settings.NSample*blksize)
    % % % %         disp('Not able to read the specified number of samples, exiting...')
    % % % %         %e_i=0.0;e_q=0.0;p_i=0.0;p_q=0.0;l_i=0.0;l_q=0.0;
    % % % %         fclose(file);
    % % % %         return
    % % % %     end
    %define index into early code vector
    tcode=(remcodephase-earlylate):codephasestep:((blksize-1)*codephasestep+remcodephase-earlylate);
    tcode2=ceil(tcode)+1;
    earlycode=ca(tcode2);
    %define index into late code vector
    tcode=(remcodephase+earlylate):codephasestep:((blksize-1)*codephasestep+remcodephase+earlylate);
    tcode2=ceil(tcode)+1;
    latecode=ca(tcode2);
    %define index into prompt code vector
    tcode=remcodephase:codephasestep:((blksize-1)*codephasestep+remcodephase);
    tcode2=ceil(tcode)+1;
    promptcode=ca(tcode2);
    
    %now compute the remainder for next time around
    %     remcodephase = (tcode(blksize) + codephasestep) - 1023.0;
    remcodephase = (tcode(blksize) + codephasestep) - numchips;
    
    %generate the carrier frequency to mix the signal to baseband
    time=(0:blksize) ./ fs;
    %get the argument to sin/cos functions
    trigarg = ((carrfreq * 2.0 * pi) .* time) + remcarrphase;
    %compute the "leftover" on sample after the last to start there next time
    remcarrphase=rem(trigarg(blksize+1),(2 * pi));
    %finally compute the signal to mix the collected data to bandband
    carrcos=cos(trigarg(1:blksize));
    carrsin=sin(trigarg(1:blksize));
    
    %generate the six standard accumulated values
    %first mix to baseband
    %     if settings.IF==0 % Baseband signal
    %         %         iBasebandSignal = I.*carrsin - Q.*carrcos;
    %         %         qBasebandSignal = I.*carrcos + Q.*carrsin;
    %         tempdatasin = I.*carrsin - Q.*carrcos;
    %         tempdatacos = I.*carrcos + Q.*carrsin;
    %     else
    %         %         qBasebandSignal = carrCos .* rawSignal;
    %         %         iBasebandSignal = carrSin .* rawSignal;
    %         tempdatacos = carrcos .* rawdata;
    %         tempdatasin = carrsin .* rawdata;
    %     end
    
    carrExp = carrsin + sqrt(-1)*carrcos;
    
    %% Generate the six standard accumulated values ---------------------------
    % First mix to baseband
    BasebandSignal = carrExp .* rawdata;
    tempdatacos = imag(BasebandSignal);
    tempdatasin = real(BasebandSignal);
    
    
    %now get early, late, and prompt values for each
    e_i(loopcnt) = sum(earlycode .* tempdatasin);
    e_q(loopcnt) = sum(earlycode .* tempdatacos);
    p_i(loopcnt) = sum(promptcode .* tempdatasin);
    p_q(loopcnt) = sum(promptcode .* tempdatacos);
    l_i(loopcnt) = sum(latecode .* tempdatasin);
    l_q(loopcnt) = sum(latecode .* tempdatacos);
    
    
    % FLL discriminator and feed-back
    %implement carrier loop discriminator (phase detector)
    if (loopcnt>1)
        dot = p_i(loopcnt-1)*p_i(loopcnt) +  p_q(loopcnt-1)*p_q(loopcnt);
        cross = p_i(loopcnt-1)*p_q(loopcnt) -  p_i(loopcnt)*p_q(loopcnt-1);
        phase = (atan2(cross,dot));    % FLL discriminator (atan2) with recovered phase-shift
        counter = counter + 1;
        if (phase>(pi/2)),
            phase = phase-pi;
            counter = 1;
        elseif (phase<-(pi/2))
            phase = phase+pi;
            counter = 1;
        end
        
        errort_carr = phase/(delta*360);
    end
    
    vec_errort_carr(loopcnt) = errort_carr;
    %implement code loop filter and generate NCO command
    carrier_nco = oldcarrier_nco + (t2_div_t1_carr * (errort_carr - olderrort_carr)) + errort_carr*delta_div_t1_carr;
    %modify code freq based on NCO command
    carrfreq =  mycarrfreq + carrier_nco;
    % N.B! freq = d(phase)/dt;
    vec_nco_carr(loopcnt) = carrier_nco;
    vec_freq_carr(loopcnt) = carrfreq;
    olderrort_carr = errort_carr;
    oldcarrier_nco = carrier_nco;
    
    %-------------------------------
    
    % DLL discriminator and feed-back
    errort_code = 0.5*(( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) - ( l_i(loopcnt)^2 + l_q(loopcnt)^2 ))/(( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) + ( l_i(loopcnt)^2 + l_q(loopcnt)^2 ));
    %errort_code = ((e_i(loopcnt) - l_i(loopcnt))*sign(p_i(loopcnt)))/(p_i(loopcnt));
    % N.B: Normalized Noncoherent early minus late power discriminator, the gain of the overall loop is set to 1 (k=1)
    
    %implement code loop filter and generate NCO command
    code_nco = oldcode_nco + (t2_div_t1 * (errort_code - olderrort_code)) + errort_code*delta_div_t1;
    
    %modify code freq based on NCO command
    codefreq = codefreq_basis - code_nco;
    
    olderrort_code = errort_code;
    oldcode_nco = code_nco;
    vec_errort_code(loopcnt) = codefreq;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     [FLL_st,freq_est] = control_FLL_v2(carrfreq-mycarrfreq,loopcnt,5);
    [FLL_st,freq_est,Freq_sum,Old_Freq_sum] = control_FLL_v2(carrfreq-mycarrfreq,loopcnt*4,0.1,Freq_sum,Old_Freq_sum,T_int);
    %     [FLL_st,freq_est] = control_FLL_v2(carrfreq-mycarrfreq,loopcnt*20,0.01);
    
    
    %     if (FLL_st ~= 0)
    % if (loopcnt*20 ==  msec)
    if (FLL_st ~= 0 || loopcnt*(T_int*1e3) ==  msec)
        %     if (FLL_st ~= 0 || loopcnt*20 ==  msec)
        carrfreq = mycarrfreq + freq_est;
        %         close all
        %         x_ms = 1:1:loopcnt;
        step_ms = (T_int*1e3);
        x_ms = step_ms:step_ms:loopcnt*step_ms;
        %         x_ms = 20:20:loopcnt*20;
        
        if plotme==1
        figure(200),subplot(121),plot(x_ms,p_i)
        grid, hold on;
        xlabel('milliseconds')
        ylabel('inphase amplitude')
        title('Prompt Inphase Correlation Results')
        subplot(122),plot(x_ms,p_q)
        grid
        xlabel('milliseconds')
        ylabel('quadrature amplitude')
        title('Prompt Quadrature Correlation Results')
        
        for i=1:length(vec_freq_carr)
            
            if (i==1),
                vec_freq_mean(i) = vec_freq_carr(i);
            else
                vec_freq_mean(i) = 1/i*vec_freq_carr(i) + ((i-1)/i)* vec_freq_mean(i-1);
            end
            
        end
        
        figure(201)
        plot (x_ms,vec_freq_carr - mycarrfreq,'r'),grid on;
        hold on
        plot (x_ms,vec_freq_mean - mycarrfreq,'g.-');
        xlabel('milliseconds');
        ylabel('Frequency [Hz]')
        legend('Increment','Average')
        title('Increment of the local carrier frequency with respect to the nominal value')
        %    figure
        %    plot (vec_errort,'k'),grid on;
        %disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);
        %mean(vet_freq)
        figure(202)
        plot(x_ms,p_i .^2 + p_q .^ 2, 'g.-')
        hold on
        grid
        plot(x_ms,e_i .^2 + e_q .^ 2, 'bx-');
        plot(x_ms,l_i .^2 + l_q .^ 2, 'r+-');
        xlabel('milliseconds')
        ylabel('amplitude')
        title('Correlation Results')
        legend('prompt','early','late')
        figure(203)
        plot(x_ms,vet_freq,'k.-'), grid on;
        xlabel('milliseconds')
        ylabel('Frequency [Hz]')
        title('Local Code Frequency')
        
        mean_freq = mean(vet_freq);
        end
        
        break;
    end;
    
    loopcnt=loopcnt+1;
    
end

% % close all
%         subplot(121),plot(p_i)
%         grid
%         xlabel('milliseconds')
%         ylabel('inphase amplitude')
%         title('Prompt Inphase Correlation Results')
%         subplot(122),plot(p_q)
%         grid
%         xlabel('milliseconds')
%         ylabel('quadrature amplitude')
%         title('Prompt Quadrature Correlation Results')
%
%         for (i=1:length(vec_freq_carr))
%
%             if (i==1),
%                 vec_freq_mean(i) = vec_freq_carr(i);
%             else
%                 vec_freq_mean(i) = 1/i*vec_freq_carr(i) + ((i-1)/i)* vec_freq_mean(i-1);
%             end
%
%         end
%
%         figure
%         plot (vec_freq_carr - mycarrfreq,'r'),grid on;
%         hold on
%         plot (vec_freq_mean - mycarrfreq,'g.-');
%         xlabel('milliseconds');
%         ylabel('Frequency [Hz]')
%         title('Increment of the local carrier frequency with respect to the nominal value')
%         %    figure
%         %    plot (vec_errort,'k'),grid on;
%         %disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);
%         mean(vet_freq)
%         figure
%         plot(p_i .^2 + p_q .^ 2, 'g.-')
%         hold on
%         grid
%         plot(e_i .^2 + e_q .^ 2, 'bx-');
%         plot(l_i .^2 + l_q .^ 2, 'r+-');
%         xlabel('milliseconds')
%         ylabel('amplitude')
%         title('Correlation Results')
%         legend('prompt','early','late')
%         figure
%         plot(vet_freq,'k.-'), grid on;
%         xlabel('milliseconds')
%         ylabel('Frequency [Hz]')
%         title('Local Code Frequency')
%         mean_freq = mean(vet_freq);
%
        fprintf('\n');

end


