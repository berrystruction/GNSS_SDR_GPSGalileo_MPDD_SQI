function [CNo_Emanuela, CNo_bluebook, CNo_SNV, CNo_MM, CNo_Bea,...
    vec_freq_carr, vet_freq, vec_nco_code, tau_m,...
    correlators,...
    N,Nallignment,ufinal,Ufinal,wFinal,yFinal,ynewFinal,dfinal,Dfinal,eFinal,powersFinal,error_vecFinal] =...
    trackcarrPLL_DLL_testjitter_4ms(file_in,settings,phaseFLL,mycarrfreq,plotme)
% function [e_i,e_q,p_i,p_q,l_i,l_q,vec_freq_carr] = trackcarrPLL_DLL_testjitter_10ms(file_in,f_sampling,prn,skp_msec,phaseFLL,mycarrfreq,msec,plotme);

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
%  e_i,e_q,p_i,p_q,l_i,l_q - the 6 output accumulator values recorded once each code period corresponding
%                            to early, late, and prompt (e, l, p) and inphase and quadrature (i, q)
%
%Requires
% m-file: cacode.m (to generate the C/A code, 1 sample/chip)
% input data file named "gpsdata.bin" of signed char type data format  "gpsdata_long.bin" of signed char type data format
% "numsamp.bin" of short type data format

% copia dei parametri in input
f_sampling=settings.f_sampling;
prn=settings.PRN; skp_msec=settings.seek_sec;
msec=settings.track_time; msecthreshold=settings.mvgAvgtime;
M=settings.M; Npoints=settings.multicorr.Npoints; maxlag=settings.multicorr.maxlag; 
version=settings.LAFversion;

% global f_sampling;
% f_sampling = 16.3676e6;

%fseek(file_in,0,'bof');

%skip through that data file to start at the appropriate sample (corresponding to code phase)
% skp_msec = rem(skp_msec,20);
% fseek(file_in,(phaseFLL-1)+round((20-skp_msec)*f_sampling*1e-3),'bof');
% IMPORTANT: mycodephase: number of samples skipped during the acquisition and processed by the FLL from the beginning of file;
%            round((20-cnt_skp)*f_sampling*1e-3) is the number of samples to skip, in order to start the phase tracking at the beginning of the next data bit.

fseek(file_in,(phaseFLL-1),'bof');  % SCHAR!!!
%fseek(file_in,(phaseFLL-1)+round((20-skp_msec)*f_sampling*1e-3),'bof');  % int8!!!
% fseek(file_in,(phaseFLL-1)*16,'bof'); % DOUBLE!!


% %get a vector with the C/A code sampled 1x/chip
% ca=cacode(prn);
% %then make 1st value be the last chip and the last value the very first chip
% ca=[ca(1023) ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca(1)];
global T_int;
% Generate the local GPS L2C code
% Local code generation
[ca, Rc] = GenerateLocCode(prn,settings.GNSS_signal);
%define number of chips in a code period
%numchips=length(ca); % 1 code period
% numchips=length(ca)*10; % 10 ms!!!
% numchips=length(ca)*20; % 20 ms!!!
numchips=length(ca)*T_int*1e3; % T_int ms!!!
%then make 1st value be the last chip and the last value the very first chip
% ca=[ca(end) ca ca(1)];
% ca=[ca(end) ca ca ca ca ca ca ca ca ca ca ca(1)];
%ca=[ca(end) ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca ca(1)]; % 20 ms

% To create the local code proportional to the integration time
ca2=ca;
for iii=1:T_int*1e3-1
    ca2=[ca2 ca];
end
ca=[ca(end) ca2 ca(1)]; 
% -------------------------------------------------------------------------

% perform various initializations
%loop counter for times through the loop
loopcnt=1;
%define sampling frequency
% fs=16.3676e6;
fs = f_sampling;

% define initial code frequency and code frequency basis of NCO in chips/sec
% codefreq = 1.0230e6;
codefreq = Rc;
% NCO values initialization
% codefreq_basis = 1.023e6;
codefreq_basis = Rc;
%define number of chips in a code period
% numchips=1023*20;
% numchips=length(ca);
%define code phase (in chips) which is "left over"
remcodephase = 0.0;
%define early-late offset (in chips)
earlylate = settings.dllCorrelatorSpacing ;...0.1

%define carrier frequency which is used over whole tracking period (ok if apprx as noncoherent DLL)
carrfreq = mycarrfreq;
carrierfreq_basis = mycarrfreq;
%define how much carrier phase is "left over" at the end of each code period to reset trig argument
remcarrphase = 0.0;
%--------------------------------------------------------------------------


% CARRIER TRACKING LOOP
zeta_carr = settings.pllDampingRatio;%1/sqrt(2);      % Butterworth
Bl_carr = settings.pllNoiseBandwidth; % 5; % 15;
%wn_carr = Bl_carr/0.53;            % Natural Frequency;
k_carr = 10; %2;  % 10;               % Gain of the overall loop
% loop filter coefficients
% t1_carr = k_carr/(wn_carr*wn_carr);
% t2_carr = 2*zeta_carr/wn_carr;

% Calcolo coefficienti del filtro
[t1_carr, t2_carr] = calcLoopCoef(Bl_carr,zeta_carr,k_carr); 

% filter values initialization
t2_div_t1_carr = t2_carr/t1_carr;
olderrort_carr = 0;
oldcarrier_nco = 0;
% delta_carr = (10e-3);               % integration time
% delta_carr = (1e-3);               % integration time
% global T_int;
delta_carr = T_int;               % integration time
delta_div_t1_carr = delta_carr/t1_carr;
vec_freq_carr = [];
vec_errort_carr = [];
vec_nco_carr =[];
%% salvo le freq di codice
vec_freq_code=[];
vec_nco_code=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CODE TRACKING LOOP
%zeta = 1/sqrt(2);      % Butterworth
zeta=settings.dllDampingRatio;  
% Bl = 2; original value 
Bl=settings.dllNoiseBandwidth;
%wn=Bl/0.53;              % Natural Frequency;
k = 1;                 % Gain of the overall loop

% loop filter coefficients
%t1 = k/(wn*wn);
%t2 = 2*zeta/wn;
% Calcolo coefficienti del filtro
[t1, t2] = calcLoopCoef(Bl,zeta,k);

% filter values initialization
t2_div_t1 = t2/t1;
olderrort_code = 0;
oldcode_nco = 0;
% delta = (10e-3);               % integration time
% delta = (1e-3);               % integration time
delta = T_int;               % integration time
delta_div_t1 = delta/t1;
blksize_vect = [];
Test_vect = [];
Discrim_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C/No estimator
global CNo_WINDOW;
global CNo_Emanuela;
CNo_Emanuela = [];
global CNo_bluebook;
CNo_bluebook = [];
global CNo_SNV;
CNo_SNV = [];
global CNo_MM;
CNo_MM = [];
global CNo_Bea;
CNo_Bea = [];

%%%% LAF initialization
alignment=true;
if alignment==true
    Nallignment=2*Npoints+1-floor(Npoints/2);
    %Nallignment=2*Npoints+1;
    N=2*Npoints+1;
else
    Nallignment=2*Npoints+1;
    N=Nallignment;
end

tott=0;%M/2;
%version=2; % if version==2, LAF standard function (covariance), version~=2
%autocorrelation method of LAF or others
if version==2
    correctNumber=Nallignment-M+1;% corretta dimensione POST LAF    
else
    correctNumber=Nallignment+M-1;% N+M-1;  autocorrelation method
    %correctNumber=Nallignment-M+1+floor(tott); % CLS-LAF
end

u1=zeros(N,1); %u1(ceil(N/2))=1;
d1=zeros(N,1); %d1(ceil(N/2))=1;

step=maxlag/Npoints; % parametri per il multicorrelatore
minlag=step;

correctDimension=floor(msec/msecthreshold/(T_int*1e3)); % dimensione della matrice (temporalmente) corretta

powersFinal=zeros(correctDimension,5);  % Powers matrix
yFinal=zeros(correctNumber,correctDimension);
ynewFinal=zeros(correctNumber,correctDimension);
wFinal=zeros(M,correctDimension);         % weights of ideal correlation
error_vecFinal=zeros(M,correctDimension);   % error vector for covariance calculation
eFinal=zeros(correctNumber,correctDimension);
Ufinal=zeros(correctNumber,M,correctDimension);           % Correlation Matrix
Dfinal=zeros(correctNumber,correctDimension);
ufinal=zeros(Nallignment,correctDimension);
dfinal=ufinal;

mvgAvgindex=1;
%%%%

tau_m=zeros(Nallignment,correctDimension);
MulticorrTs=false; % per abilitare/disabilitare il multicorrelatore a passi di Ts
%start a timer
tic;
% while less than the number of specified code periods, continue processing
% while (loopcnt <=  msec/10)
% while (loopcnt <=  msec/20)
while (loopcnt*(T_int*1e3) <=  msec)
    
    %since it can be time consuming, do a periodic display of current stage
    %also a good point to display debugging information
    %    if (rem(loopcnt,10)==0)
    if (rem(loopcnt*(T_int*1e3),100)==0)
        %        disp(['   Currently on interation:',int2str(loopcnt*10),' of:',int2str(msec)])
        disp(['   Currently Instant:',int2str((loopcnt*(T_int)+skp_msec)*1e3),' of:',int2str(skp_msec*1e3+msec) ' ms'])
    end
    
    %    ? update on the basis of the code loop
    %    %read in the blksize from the previous generated code tracking
    %    blksize = fread(fid2,1,'short');
    
    %read in the appropriate number of samples of data and prompt prn code to process this interation
    %    [rawdata,scount] = fread(file_in,blksize,'schar');
    %    rawdata=rawdata';
    %
    %update the phasestep based on code freq (variable) and sampling frequency (fixed)
    codephasestep=codefreq/fs;
    vet_freq(loopcnt) = codefreq;
    
    %find the size of a "block" or code period in whole samples
    %note - not starting from zero, but rather where you left off from last code period
    %this could be done using the CA code generator with sampling freq, code freq, and offset inputs
    %but this will be faster
    blksize=ceil((numchips-remcodephase)/codephasestep);
    blksize_vect(loopcnt) = blksize;
    %read in the appropriate number of samples to process this interation
    [rawdata,scount] = fread(file_in,blksize,settings.dataType);
    %    [rawdata,scount] = fread(file_in,blksize,'double'); % DOUBLE!!
    rawdata=rawdata';
    %if did not read in enough samples, then could be out of data - better exit
    if (scount ~= blksize)
        disp('Not able to read the specified number of samples, exiting...')
        % e_i=0.0;e_q=0.0;p_i=0.0;p_q=0.0;l_i=0.0;l_q=0.0; % provo a
        % lasciarla commentata e salvo i 5s di correlatori
        fraction=size(e_i,2);
        correlators=struct('ei',e_i(1:fraction),'eq',e_q(1:fraction),'pi',p_i(1:fraction),'pq',p_q(1:fraction),'li',l_i(1:fraction),'lq',l_q(1:fraction));
        %%% vettore correzione di frequenza di codice
        %vec_nco_code=[0 vec_nco_code];
        
        %% ritorno le variabili con la dimensione corretta proporzionata al tempo di elaborazione
        correctDimension=floor((loopcnt*(T_int*1e3))/msecthreshold);
        powersFinal=powersFinal(1:correctDimension,:);  % Powers matrix
        yFinal=yFinal(:,1:correctDimension);
        ynewFinal=ynewFinal(:,1:correctDimension);
        wFinal=wFinal(:,1:correctDimension);         % weights of ideal correlation
        error_vecFinal=error_vecFinal(:,1:correctDimension);   % error vector for covariance calculation
        eFinal=eFinal(:,1:correctDimension);
        Ufinal=Ufinal(:,:,1:correctDimension);           % Correlation Matrix
        Dfinal=Dfinal(:,1:correctDimension);
        %fclose(file_in);
        return
    end
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
    remcodephase_old = remcodephase;
    
    %    remcodephase = (tcode(blksize) + codephasestep) - 1023.0*20;
    remcodephase = (tcode(blksize) + codephasestep) - numchips;
    
    Test_vect(loopcnt) = (blksize-1)*(1/fs)  + (codephasestep-remcodephase)*(1/codefreq) + remcodephase_old*(1/codefreq);
    %Test_vect(loopcnt) = (blksize-1)*(1/fs)  - remcodephase*(1/codefreq);
    
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
    tempdatacos = carrcos .* rawdata;
    tempdatasin = carrsin .* rawdata;
    
    %now get early, late, and prompt values for each
    e_i(loopcnt) = sum(earlycode .* tempdatasin);
    e_q(loopcnt) = sum(earlycode .* tempdatacos);
    p_i(loopcnt) = sum(promptcode .* tempdatasin);
    p_q(loopcnt) = sum(promptcode .* tempdatacos);
    l_i(loopcnt) = sum(latecode .* tempdatasin);
    l_q(loopcnt) = sum(latecode .* tempdatacos);
    
    
    %% -------------------------------
    
    % PLL discriminator and feed-back
    %implement carrier loop discriminator (phase detector)
    errort_carr = atan ( p_q(loopcnt)/p_i(loopcnt))/ (2.0 * pi);    % discriminator as arc-tan
    
    vec_errort_carr(loopcnt) = errort_carr;
    %implement code loop filter and generate NCO command
    carrier_nco = oldcarrier_nco + (t2_div_t1_carr * (errort_carr - olderrort_carr)) + errort_carr*delta_div_t1_carr;
    %modify code freq based on NCO command
    %delta_freq = code_nco - oldcode_nco;
    carrfreq =  mycarrfreq + carrier_nco;
    % N.B! freq = d(phase)/dt;
    vec_nco_carr(loopcnt) = carrier_nco;
    vec_freq_carr(loopcnt) = carrfreq;
    olderrort_carr = errort_carr;
    oldcarrier_nco = carrier_nco;
    
    %% -------------------------------
    
    % DLL discriminator and feed-back
    errort_code =0.5*((( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) - ( l_i(loopcnt)^2 + l_q(loopcnt)^2 ))/(( e_i(loopcnt)^2 + e_q(loopcnt)^2 ) + (l_i(loopcnt)^2 + l_q(loopcnt)^2 )));
    %errort_code = ((e_i(loopcnt) - l_i(loopcnt))*sign(p_i(loopcnt)))/(p_i(loopcnt));
    % N.B: Normalized discriminator, the gain of the overall loop is set to 1 (k=1)
    Discrim_out(loopcnt) = errort_code;
    
    %implement code loop filter and generate NCO command
    code_nco =  (t2_div_t1 * (errort_code - olderrort_code)) + errort_code*delta_div_t1 + oldcode_nco;
    %code_nco =  errort_code - olderrort_code + errort_code*0 +
    %oldcode_nco; perde aggancio

    %modify code freq based on NCO command
    codefreq = codefreq_basis - code_nco;
    % salvo i valori
    vec_nco_code(loopcnt)=code_nco;
    %vec_freq_code=codefreq;
        
    
    %% Multicorrelator (2 versions available)
    if MulticorrTs==false
        [new_u, new_d]=Multicorrelator(ca,tempdatacos,tempdatasin,codephasestep,remcodephase_old,...
            blksize,minlag,step,maxlag);
    else
        % code matrix
        offsetXcorrpoints=0; % deve essere un valore pari (al momento)
        % code_m =codematrixgen(promptcode,Npoints,offsetXcorrpoints);
        [new_u, new_d]=MulticorrelatorTs(promptcode,tempdatacos,tempdatasin,Npoints,...
            offsetXcorrpoints,1/f_sampling,remcodephase_old,codefreq);
    end
    tau=[-(Nallignment-1)/2:(Nallignment-1)/2]*step/codefreq;
    %tau=[-(length(new_d)-1)/2:(length(new_d)-1)/2]'*step/codefreq;
    tau_m(:,loopcnt)=tau;
    %%%
    if rem(loopcnt,msecthreshold)==0
        
        % add new correlation calculation  
        u1=(u1+new_u)/msecthreshold;
        d1=(d1+new_d)/msecthreshold;
        
        %% Normalization
        %         u=u/max(u);
        %         d=d/max(d);
        
        %% Alignment
        if alignment==true
            [~, indexmax]=max(u1);
            u=MaxAlignment(u1,Nallignment,indexmax);
            %[~, indexmax]=max(d1);
            d=MaxAlignment(d1,Nallignment,indexmax);
        else
            u=u1;
            d=d1;
        end

        
        %% LAF procedure
        [LAF]=LAF_SQM(u,d,M,version,tott);%remcodephase/codefreq)
        %% save return parameters
        powersFinal(mvgAvgindex,:)=LAF.powers(1,:);  % Powers matrix
        yFinal(:,mvgAvgindex)=LAF.y;
        ynewFinal(:,mvgAvgindex)=LAF.ynew;
        wFinal(:,mvgAvgindex)=LAF.w;         % weights of ideal correlation
        error_vecFinal(:,mvgAvgindex)=LAF.error_vec;   % error vector for covariance calculation
        eFinal(:,mvgAvgindex)=LAF.e;
        Ufinal(:,:,mvgAvgindex)=LAF.U;
        Dfinal(:,mvgAvgindex)=LAF.D;
        ufinal(:,mvgAvgindex)=u;
        dfinal(:,mvgAvgindex)=d;
        mvgAvgindex=mvgAvgindex+1;
        
        %% azzeramento ogni mvgAvg time?
        u1=zeros(N,1); %u(ceil(N/2))=1;
        d1=zeros(N,1); %d(ceil(N/2))=1;
        
    else  
        u1=u1+new_u;
        d1=d1+new_d;     
    end
    %%
    
    olderrort_code = errort_code;
    oldcode_nco = code_nco;
    vec_errort_code(loopcnt) = errort_code;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     fwrite(fid_1,p_i(loopcnt),'double');
    %     fwrite(fid_2,p_q(loopcnt),'double');
    
    PLL_st = control_PLL_v2_EF(p_i(loopcnt),p_q(loopcnt),loopcnt,35,f_sampling);
    
    if (PLL_st ~= 0)
        disp('PLL not locked, C/No below the threshold, the channel might be reassigned to another satellite, or try with a new fast acquisition...')
        break;
    end;
    
    %increment loop counter
    loopcnt=loopcnt+1;
    
end
% fclose(file_in);

fraction=size(e_i,2);
correlators=struct('ei',e_i(1:fraction),'eq',e_q(1:fraction),'pi',p_i(1:fraction),'pq',p_q(1:fraction),'li',l_i(1:fraction),'lq',l_q(1:fraction));
%%% vettore correzione di frequenza di codice
%vec_nco_code=[0 vec_nco_code];

%%
disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec'])

if (plotme == 1)
    %    close all
    
    %     x_axe_second = 0.01:0.01:length(p_i)/100;
    step_ms = (T_int*1e3);
    x_axe_second = step_ms:step_ms:length(p_i)*step_ms;
    
    subplot(121),plot(x_axe_second,p_i,'o-')
    grid
    xlabel('milliseconds')
    ylabel('In-phase amplitude')
    title('Prompt In-phase Correlation Results')
    subplot(122),plot(x_axe_second,p_q,'o-')
    grid
    xlabel('milliseconds')
    ylabel('Quadrature amplitude')
    title('Prompt Quadrature Correlation Results')
    
    figure
    plot(p_i,p_q,'g.'),title('Prompt I & Q'), xlabel('I'), ylabel('Q'),grid on;
    
    %
    %    for (i=1:length(vec_freq_carr))
    %         if (i<=500),
    %             vect_local = vec_freq_carr(1:i);
    %             vec_freq_mean(i) = sum(vect_local)/length(vect_local);
    %         else
    %             vect_local = vec_freq_carr(i-499:i);
    %             vec_freq_mean(i) = sum(vect_local)/length(vect_local);
    %         end
    %    end
    for i=1:length(vec_freq_carr)
        if (i<=25),
            vect_local = vec_freq_carr(1:i);
            vec_freq_mean(i) = sum(vect_local)/length(vect_local);
        else
            vect_local = vec_freq_carr(i-24:i);
            vec_freq_mean(i) = sum(vect_local)/length(vect_local);
        end
    end
    figure
    plot (x_axe_second,vec_freq_carr - mycarrfreq,'r'),grid on;
    hold on
    plot (x_axe_second,vec_freq_mean - mycarrfreq,'g.-');
    xlabel('milliseconds');
    ylabel('Frequency [Hz]')
    title('Increment of the local carrier frequency with respect to the nominal value')
    %    figure
    %    plot (vec_errort,'k'),grid on;
    % disp(['Time required to process ',int2str(msec),' msec of data is:',num2str(toc),'sec']);
    % mean(vet_freq)
    
    figure(100)
    subplot(212),plot(x_axe_second, (p_i .^2 + p_q .^ 2), 'g.-'),grid on;
    hold on
    subplot(212),plot(x_axe_second, (e_i .^2 + e_q .^ 2), 'bx-');
    subplot(212),plot(x_axe_second, (l_i .^2 + l_q .^ 2), 'r+-');
    grid on
    xlabel('milliseconds');
    ylabel('amplitude');
    title('Correlation Results');
    legend('prompt','early','late');
    
    %     x_axe_second_CNo = 0.5:0.5:length(CNo_bluebook)/2;
    x_axe_second_CNo = (0:(length(CNo_MM)-1))./(1./CNo_WINDOW)+500;
    subplot(211), plot(x_axe_second_CNo,CNo_MM,'k.-'),grid on;
    ylabel('C/No[dBHz]');
    xlabel('Time [ms]')
    
    figure(101),
    subplot(211),plot(x_axe_second,vec_errort_carr,'k*-'),hold on,grid on
    title('Phase residual error');
    subplot(212),plot(x_axe_second,(e_i-l_i)./(2*p_i),'k-'),hold on,grid on
    title('Pseudorange tracking error [chip]');
    
    figure
    plot(x_axe_second,vet_freq,'k.-'), grid on;
    xlabel('milliseconds')
    ylabel('Frequency [Hz]')
    title('Local Code Frequency')
    mean_freq = mean(vet_freq);
    
    
    figure
    plot(x_axe_second,vec_errort_carr,'r*-'),hold on,grid on
    title('Phase residual error')
    
    figure
    plot(x_axe_second,vec_errort_code,'m*-'),hold on,grid on
    title('Code residual error')
    
    global nominalfreq;
    figure
    plot (x_axe_second,vec_freq_carr - nominalfreq,'g.-');hold on,grid on;
    xlabel('milliseconds');
    ylabel('Doppler Frequency [Hz]')
    title('Estimated Doppler frequency')
    
end



end

