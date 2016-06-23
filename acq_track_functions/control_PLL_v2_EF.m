
% Control PLL.
% This program control the PLL status and estimates the C/No ratio.


function [PLL_status,FIFO_IP_aus,FIFO_QP_aus,CN0_aus] = control_PLL_v2_EF(act_IP,act_QP,loop,thres,f_sampling,T_int,CNo_WINDOW,FIFO_IP,FIFO_QP)

samplesPDI = ceil(f_sampling*T_int);


%%%%%EF Beq = 2*(f_sampling/samplesPDI);  % real signal component + complex noise estimation
Beq = (f_sampling/samplesPDI);  % real signal component + complex noise estimation

% Update the FIFO structure used to estimate the C/No

FIFO_IP_aus = [act_IP, FIFO_IP(1:end-1)];
FIFO_QP_aus = [act_QP, FIFO_QP(1:end-1)];

PLL_status = 0;
%index = 0;

% if (rem(loop*10,1000)==0) % Estimate the C/No EVERY 1 s
% if (rem(loop*10,500)==0) % Estimate the C/No EVERY 0.5 s
%
% if ((loop*10>1500)&&(rem(loop*10,CNo_WINDOW)==0)) % Skip 1.5 s and Estimate the C/No EVERY CNo_WINDOW ms
% if ((loop*10>1500)&&(rem(loop*10,1000)==0)) % Skip 1.5 s and Estimate the C/No EVERY 1 s
% if ((loop*10>1500)&&(rem(loop*10,500)==0)) % Skip 1.5 s and Estimate the C/No EVERY 500 ms
% if ((loop*10>1500)&&(rem(loop*10,100)==0)) % Skip 1.5 s and Estimate the C/No EVERY 100 ms
%global T_int;
%if ((loop*(T_int*1e3)>500)&&(rem(loop*(T_int*1e3),CNo_WINDOW)==0)) % Skip .5 s and Estimate the C/No EVERY CNo_WINDOW ms
if (rem(loop*(T_int*1e3),CNo_WINDOW)==0) % Estimate the C/No EVERY CNo_WINDOW ms
    
    %index = index + 1;
    %*** Estimate the C/No ratio - Emanuela's Method ***
    %disp('Here estimate the C/No and check the lock of the tracking loops')
    
    Pn = 2*var(FIFO_QP_aus);
    Ptot = var(FIFO_IP_aus + 1i*FIFO_QP_aus);
    Ps = Ptot-Pn;
    if Ps > 0
        SNR_RSCN = Ps/Pn;
    else
        SNR_RSCN = 1;     % in case of Pn<0, it means that a phase rotation occured.  The estimator is not able to provide a correct value of SNR,
        % which is set equal to 1 in order to have 0 in the C/No.
    end
    
    CNo = 10*log10(SNR_RSCN * Beq);
    CN0.CNo_Emanuela = CNo;%[CN0.CNo_Emanuela, CNo];
    % ***************************************************
    %*** Estimate the C/No ratio - Blue Book Method ***
    %disp('Here estimate the C/No and check the lock of the tracking loops')
    %     MSamplesPerBit = 2; % GPS
    MSamplesPerBit = 1; % Galileo E1b 4 ms
    
    K = length(FIFO_IP_aus)/MSamplesPerBit; % represents the number of bits used in the estimation
    HTotSampleNumber = K*MSamplesPerBit;
    
    I_sam = FIFO_IP_aus(1:HTotSampleNumber);
    Q_sam = FIFO_QP_aus(1:HTotSampleNumber);
    
    CMatrix = reshape(I_sam,MSamplesPerBit,K) + 1i*reshape(Q_sam,MSamplesPerBit,K);
    
    WBP = sum(abs(CMatrix).^2);
    
    NBP = (sum(real(CMatrix))).^2 + (sum(imag(CMatrix))).^2;
    
    NP = NBP./WBP;
    
    mNP = mean(NP);
    
    CNo_PRM = (mNP-1)/(MSamplesPerBit-mNP)/10e-3;
    
    if CNo_PRM < 0
        CNo_PRM = NaN;
    end
    
    CNo = 10*log10(CNo_PRM);
    CN0.CNo_bluebook =CNo;% [CN0.CNo_bluebook, CNo];
    % **************************************************
    %*** Estimate the C/No ratio - SNV Method ***
    
    %%%%%EF Ps = mean(abs(FIFO_IP))^2 + mean(abs(FIFO_QP))^2;
    Ps = mean(abs(FIFO_IP_aus))^2;
    Ptot = mean(abs(FIFO_IP_aus+1i*FIFO_QP_aus).^2);
    Pn = Ptot-Ps;
    
    if Pn > 0
        SNR_SNV = Ps/Pn;
    else
        SNR_SNV = NaN;
    end
    
    CNo = 10*log10(SNR_SNV * Beq); %;
    %disp([ 'CN0= ' num2str(CNo) 'dBHz']);
    CN0.CNo_SNV = CNo;%[CN0.CNo_SNV, CNo];
    
    % ***************************************************
    %*** Estimate the C/No ratio - Moments Method ***
    
    m2 = mean(abs(FIFO_IP_aus+1i*FIFO_QP_aus).^2);
    m4 = mean(abs(FIFO_IP_aus+1i*FIFO_QP_aus).^4);
    Ps = sqrt(2*m2^2 - m4);
    Pn = m2-Ps;
    if Pn > 0
        SNR_MM = Ps/Pn;
    else
        SNR_MM = NaN;
    end
    
    CNo = 10*log10(SNR_MM * Beq);
    CN0.CNo_MM =CNo;% [CN0.CNo_MM, CNo];
    
    % **************************************************
    %*** Estimate the C/No ratio - Beaulieu ***
    
    %%%%%EF
    
    Nsym = length(FIFO_IP_aus)-1;
    X_re = FIFO_IP_aus(1:end-1);
    Y_im = FIFO_IP_aus(2:end);
    func2 = sum((abs(X_re) - abs(Y_im)).^2 ./ (X_re.^2 + Y_im.^2));
    if func2 > 0
        SNR_Bea = Nsym/func2/2;
    else
        SNR_Bea = NaN;
    end
    CNo = 10*log10(SNR_Bea * Beq);
    CN0.CNo_Bea =CNo;% [CN0.CNo_Bea, CNo];
    % **************************************************
    
    
    CN0_aus=CN0;
    CN0_aus.CN0ready=true;
    
else
%     CN0_aus=[];
%     CN0aus.CNo_Emanuela=0;
%     CN0_aus.CNo_bluebook=0;
%     CN0_aus.CNo_SNV=0;
%     CN0_aus.CNo_MM=0;
%     CN0_aus.CNo_Bea=0;
    CN0_aus.CN0ready=false; % in roder to disable saving parameter
end

end