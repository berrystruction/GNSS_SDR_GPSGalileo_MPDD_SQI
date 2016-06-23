% FFT Acquisition in Time domain

function [doppler, code, st, C,peakMetric] = signal_acquisition_FFT(file,T_int,settings,prn)%(file,sampling_freq,IF,prn,T_int,seek_sec)

sampling_freq=settings.f_sampling;
IF=settings.nominalfreq;
seek_sec=settings.seek_sec;

%wbh = waitbar(0,'Signal acquisition. Please wait...');
j=sqrt(-1);
% T_int = 4e-3; %1e-3; %20e-3;     % [s] coherent integration time
% T_int = CodeLen/Rc;    % [s] coherent integration time = 1 primary code period (= 4 ms)
%global T_int;
Non_Coh_Sums = 100; %3; %5; %11; %20 %50 %100; % Non coherent summations (>= 1)

Dopplerstep  = settings.acqSearchStep; %250; %25; %50; %100; %500; % [Hz]
DopplerRange = settings.acqSearchBand*1e3; %8000; %1000; %2000; %200; %5000; %2000; % [Hz]

FD_vect= -DopplerRange:Dopplerstep:DopplerRange;
Doppler_CENTER = 0; %1000 %0 %2000
FD_vect = FD_vect + Doppler_CENTER; %%%%%%%%%%%%%%%%%%%%%%%%%

% FD_vect= -20000:Dopplerstep:20000;
% FD_vect= -5000:Dopplerstep:5000;
% FD_vect= -2000:Dopplerstep:2000;
% FD_vect= -1000:Dopplerstep:1000;

acq_metric = settings.acqThreshold; % Acquisition metric
% acq_metric = 1.1; % Acquisition metric


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local code generation
[Loc, Rc] = GenerateLocCode(prn,settings.GNSS_signal);

CodeLen=length(Loc);

% T_int = CodeLen/Rc;    % [s] coherent integration time = 1 primary code period (= 4 ms)

% Loc = [Loc Loc(1)];  % add the first chip at the end of the code vector

num_samples = sampling_freq/Rc;         % Number of samples per chip.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


status = fseek(file,ceil(sampling_freq*seek_sec),-1);
if (status)
    error('fseek error')
end

%fseek(file,ceil(sampling_freq*0.1),-1);
% [data,scount] = fread(file,4194304,'int8'); %...solitamente leggere schar, da NordNav
% data = data';

% Nsamples = ceil(T_int*sampling_freq+1)*(Non_Coh_Sums+1); % samples to be read
Nsamples = floor(T_int*sampling_freq)*(Non_Coh_Sums); % samples to be read
% [data,scount] = fread(file,Nsamples,'uchar'); % INPUT SAMPLES READING
% [data,scount] = fread(file,Nsamples,'double'); % INPUT SAMPLES READING
[data,scount] = fread(file,Nsamples,'int8'); % INPUT SAMPLES READING
if (scount~=Nsamples)
    disp('Not able to read the specified number of samples, exiting...')
    %fclose(file);
    return
end

% data = data'-127; 
data = data'; 

%Ts = 1/sampling_freq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=1; % Downsampling

Fs=sampling_freq/K;

fif=rem(IF,Fs/K);

data=data(1:K:end);

% N=floor(Fs*CodeLen/Rc)+1;
% N=floor(Fs*CodeLen/Rc)
N=floor(Fs*T_int);       % Number of samples for each coherent integration

N_tot = N * (Non_Coh_Sums); % Total number of processed samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFT acquisition

idx=1;

% Ricampiono il codice locale (NB: il Doppler sul codice è trascurato)
n=0:N_tot-1;
ind_cod = mod(floor(n*Rc/Fs),CodeLen)+1;
SigLOC_tot=Loc(ind_cod);

C = zeros(length(FD_vect),N);

for ind_FD= 1:length(FD_vect)  
        
    FD = FD_vect(ind_FD);
    
    corr=zeros(1,N)+j*zeros(1,N);
          
%     % Ricampiono
%     
    k=0:N-1;
%     SigLOC=Loc(floor(k*Rc/Fs)+1);
%     
%     % FFT codice locale e complesso coniugato
%     SigLOCFFT=conj(fft(SigLOC,N));
    
    argx=2*pi*(fif+FD)/Fs;
    carrI=cos(argx*k);
    carrQ=sin(argx*k);
    
%     for M=0:5
    for M=0:(Non_Coh_Sums-1)
        
        % FFT codice locale e complesso coniugato
        SigLOC=SigLOC_tot(N*M+1:N*M+N);
        SigLOCFFT=conj(fft(SigLOC,N));
        
        % Dati ricevuti
        SigIN=data(N*M+1:N*M+N);
        
        % Demodulo
        
        I=SigIN.*carrI;
        Q=SigIN.*carrQ;
        
        SigINIQ=I+j*Q;
        
        corr=corr+abs(ifft(fft(SigINIQ,N).*(SigLOCFFT)));
        
    end
    
    C(idx,:)=corr;
    idx=idx+1;

    %waitbar(ind_FD/length(FD_vect),wbh)
    
end

%wbh.delete;

% --- Find the main peak in the correlation floor and the corresponding frequency bin index
[bb, ind_mixf] = max(max(C,[],2),[],1);
[bb, ind_mixc] = max(max(C));

if (ind_mixc < ceil(num_samples/K)),
    vect_search_peak2 = [zeros(1,2*ceil(num_samples/K)), C(ind_mixf,(2*ceil(num_samples/K)):end)];    
elseif (ind_mixc > (length(C(ind_mixf,:)-2*ceil(num_samples/K))))
    vect_search_peak2 = [C(ind_mixf,1:(end-2*ceil(num_samples/K)):end), zeros(1,2*ceil(num_samples/K))];  
else
    vect_search_peak2 = [C(ind_mixf,1:(ind_mixc-ceil(num_samples/K))),zeros(1,2*ceil(num_samples/K)-1),C(ind_mixf,(ind_mixc+ceil(num_samples/K)):end)];
end

% --- Find the second highest peak in the correlation floor
second_peak = max(vect_search_peak2);


% --- compare the acquisition metric to a predefined threshold

peakMetric=bb/second_peak;

if ((bb/second_peak) > acq_metric)
    fprintf('...acquired satellite\n')
    code = ceil(sampling_freq*seek_sec)+((ind_mixc-1)*K);% - 20;
    % 20 is a sistematic error, which simulates the number of offset samples, which makes the acquisition estimation not correct, since it took time (not negligible).
    % The effect on the carrier frequency can be neglected, no problems for that!
      
    doppler = (IF+FD_vect(1)) + (ind_mixf-1)*Dopplerstep;
       
    st = 1;
else
    fprintf('...no satellite\n')
    code = 0;
    doppler = 0;
    st = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
figure,
plot([FD_vect],C(:,ind_mixc))
xlabel('Frequency Bins [Hz]')
ylabel('Value')
title('Peak varying Doppler bin')
grid on

figure;
plot(k,C(ind_mixf,:))
xlabel('Code Delay Bins [samples]')
ylabel('Value')
title(' Peak varying Code phase')
grid on

if numel(C)<1e6
    figure
    surf(k,FD_vect,C/max(max(C)))
    shading interp
    axis([k(1) k(end) FD_vect(1) FD_vect(end) 0 1])
    title('Normalized Search Space')
    xlabel('Code Delay Bins [samples]')
    ylabel('Frequency Bins [Hz]')
    zlabel('Normalized correlation')
end

%close(wbh)
end