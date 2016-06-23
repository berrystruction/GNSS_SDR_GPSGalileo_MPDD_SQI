% SERIAL Acquisition in Time domain (acquisition refinement)

function [doppler, code] = signal_acquisition_SERIAL(file,fs,prn,code_Ph,mix_freq)

% file: is the original file, where the raw samples of the signal are stored;
% fs: sampling frequency;
% prn: PRN to track;
% code_Ph: code_phase in samples, where the main peak has been detected in the acquisition phase;
% mix_freq: mixing frequency where the main peak has been detected in the acquisition phase

wbh = waitbar(0,'Serial acquisition. Please wait...');

T_int = 1e-3; %20e-3;     % [s] coherent integration time
Non_Coh_Sums = 5; %11; %11; %50; %20; %100; % Non coherent summations (>= 1)

Chips_Pull = 3; %1; %2;   % Number of chips where the Pull_in phase must search for correlation peaks

Dopplerstep  = 100; %1; %25; %50; %100; %500; % [Hz]
DopplerRange = 2000; %25; %50; %200; %5000; %2000; % [Hz]

FD_vect= -DopplerRange:Dopplerstep:DopplerRange;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local code generation
[Loc, Rc] = GenerateLocCode(prn);

CodeLen=length(Loc);

% Loc = [Loc Loc(1)];  % add the first chip at the end of the code vector

num_samples = fs/Rc;         % Number of samples per chip.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File reading
Nsamples = floor(T_int*fs)*(Non_Coh_Sums) + code_Ph + (ceil(Chips_Pull*num_samples)); % samples to be read
fseek(file,0,-1);
% [data,scount] = fread(file,Nsamples,'uchar'); % INPUT SAMPLES READING
[data,scount] = fread(file,Nsamples,'int8'); % INPUT SAMPLES READING
%[data,scount] = fread(file,Nsamples,'double'); % INPUT SAMPLES READING
if (scount~=Nsamples)
    disp('Not able to read the specified number of samples, exiting...')
    fclose(file);
    return
end

% data = data'-127; 
data = data'; 

data = data(code_Ph - (floor(Chips_Pull/2*num_samples)):end);

bb = length(data);               % length of the new input vector, after the skip of the samples

Ts = 1/fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=1; % Downsampling
% 
Fs=fs/K;
% 
% fif=rem(IF,Fs/K);
% 
% data=data(1:K:end);

N=floor(Fs*T_int);       % Number of samples for each coherent integration

N_tot = N * (Non_Coh_Sums); % Total number of processed samples



% % %------------------------------------------------------------------
% % Code_rate_L5 = 10.23*1e6;              % Nominal GPS L5 code rate.
% % CodeLenL5 = 10230;                     % L5 primary Code Lenght 
% % num_samples = fs/Code_rate_L5;         % Number of samples per chip.
% % 
% % fseek(file,0,-1);
% % [gpsdata,scount] = fread(file,500000,'int8'); %...solitamente leggere schar, da NordNav
% % gpsdata = gpsdata(code_Ph - ((Chips_Pull/2)*ceil(num_samples)):end);
% % 
% % bb = length(gpsdata);               % length of the new input vector, after the skip of the samples
% % 
% % % Generate and sample the local L5 primary code 
% % % generate the GPS L5 primary code - real takes the I part, imag takes the Q part 
% % Loc_I = real(GenModulation_GPS_L5BU(prn));  
% % Loc_I = [Loc_I, Loc_I(1)]; % add the first chip at the end of 10230 chips.
% % 
% % % Number of samples per 1 code period
% % N = floor(fs*CodeLenL5/Code_rate_L5)+1;
% % k = 0:N-1;
% % SigLoc_I = Loc_I(floor(k*Code_rate_L5/fs)+1);
% % %SigLoc_I_half = SigLoc_I(1:floor(length(SigLoc_I)/2)); % do partial correlation, take 1/2 ms of local code


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Serial acquisition

% Ricampiono il codice locale (NB: il Doppler sul codice è trascurato)
n=0:N_tot-1;
ind_cod = mod(floor(n*Rc/Fs),CodeLen)+1;
SigLOC_tot=Loc(ind_cod);

corr = zeros(length(FD_vect),floor(num_samples*Chips_Pull));

% for FD=-1500:250:1500,
for ind_FD= 1:length(FD_vect)  
        
    FD = FD_vect(ind_FD);

    indb = mix_freq + FD;

    timev=0:1/Fs:(bb-1)/Fs;  %create a sample time vector

    sinv = sin(2*pi*indb*timev);  %create a sin vector
    cosv = cos(2*pi*indb*timev);  %create a cos vector

    datasin = sinv .* data;  %mix data to baseband with sin vector
    datacos = cosv .* data;  %mix data to baseband with cos vector

    for indc=0:((floor(num_samples*Chips_Pull/2))*2-1)  % correlation shift done every sample
        corr_tmp = 0;
        for M=0:(Non_Coh_Sums-1) % Non coherent accumulations
            sinresult = sum(datasin(1+indc+N*M:N+N*M+indc) .* SigLOC_tot(N*M+1:N*M+N));  %correlation with sin downconverted
            cosresult = sum(datacos(1+indc+N*M:N+N*M+indc) .* SigLOC_tot(N*M+1:N*M+N));  %correlation with cos downconverted
            corr_tmp = corr_tmp + sqrt( sinresult * sinresult + cosresult * cosresult );  %remember trig identify: i^2+q^2 to get mag
        end
        corr(ind_FD,indc+1) = corr_tmp;
    end
    
    waitbar(ind_FD/length(FD_vect),wbh)
end

vet_maxs = max(corr);
bb = max(vet_maxs);
[ind_mixf, ind_code] = find(corr == bb);

% --- return the correct code phase and mixing frequency
code = code_Ph - (floor(Chips_Pull/2*num_samples)) + (ind_code-1);

% doppler = mix_freq + (ind_mixf-1)*100-1500;
doppler = (mix_freq+FD_vect(1)) + (ind_mixf-1)*Dopplerstep;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
figure
k = (0:(floor(num_samples*Chips_Pull)-1)) - (floor(Chips_Pull/2*num_samples));
surf(k,FD_vect,corr/max(max(corr)))
shading interp
axis([k(1) k(end) FD_vect(1) FD_vect(end) 0 1])
title('Normalized Search Space (acquisition refinement)')
xlabel('Code Delay Bins [samples]')
ylabel('Frequency Bins [Hz]')
zlabel('Normalized correlation')

close(wbh)

return

