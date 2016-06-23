function [all_is_ok, PRNindexforLocalTest,TGlobalTest,TlocalTest,Globalthresh,PL]=FaultDetection(xyzdt,RAIMparam,settings,cov_noise,SQI,KFsolution)
% Fault detection and exclusion function

% CN0sorted=[50 46 46 45 43 43 44 45 42];
% % Uso anche la formula di Lachapelle (for ligthly degraded signal condition)
% a=10;  % [m^2]
% b=150^2; % [m^2*Hz]
% sigma_i=a+b*10.^(-CN0sorted/10);
Nmeasurement=size(RAIMparam.w,1);
nu=Nmeasurement-4;

%for ausindex=1:Nmeasurement, cov_noise(ausindex,ausindex)=cov_noise(ausindex,ausindex)*SQI(ausindex)/sum(SQI); end

W=(cov_noise)^(-1);
%W=1/cov_noise(1,1)*eye(size(W));

HWH=(RAIMparam.H'*W*RAIMparam.H)^(-1);

covw=cov_noise-RAIMparam.H*HWH*RAIMparam.H';
if nargin==6
    % Test parameters for the KF
    d=RAIMparam.w;
    S=RAIMparam.H*RAIMparam.aprioriPk(1:4,1:4)*RAIMparam.H'+cov_noise;
    TGlobalTest=sqrt(d'*S^(-1)*d);
    
    wnormalized=d./sqrt(diag(covw));
else
    % Test parameters for the LS and SQI
    %RAIMparam.w=RAIMparam.w.*(1./SQI)/sum(1./SQI);
    wnormalized=RAIMparam.w./sqrt(diag(covw));
    TGlobalTest=RAIMparam.w'*W*RAIMparam.w;
end


[TlocalTest, PRNindexforLocalTest]=max(abs(wnormalized)); % observation with the largest abs value is tested

Pmd = settings.RAIM.Pmd;%0.005; %0.2
alpha0 = settings.RAIM.PfaLocal;%0.001; % 0.107
%alpha = 0.001; LA CALCOLO

lambda=(norminv(1-alpha0/2,0,1)+norminv(1-Pmd,0,1))^2;

% Globalthresh calculation
%xbins=0:.1:50;
%chi2standard=chi2pdf(xbins,nu);
%diff=abs(ncx2pdf(xbins,nu,lambda)-chi2standard);
%[minn, minpos]=min(diff(2:end)); % escludo lo zero nell'origine
%Globalthresh=xbins(minpos);

% Revised calculation for local test consideration
Globalthresh=icdf('ncx2',Pmd,nu,lambda);
% alpha, as the global Pfa
alpha=1-cdf('chi2',Globalthresh,nu);
settings.RAIM.PfaGlobal=alpha;


all_is_ok=1; % Initialization

% Global test
if TGlobalTest>Globalthresh
    GlobalTestresults=1; % failed
    % Performing local test...
    if TlocalTest<=norminv(1-alpha0/2,0,1);
        PRNindexforLocalTest=0; % solution unreliable
        all_is_ok=0;
    else
        all_is_ok=-1;
    end
    PL.HRMS=NaN;
    PL.VRMS=NaN;
    PL.PRMS=NaN;
    PL.HSlopeMax=NaN;
    PL.VSlopeMax=NaN;
    PL.HPL=NaN;
    PL.VPL=NaN;
    PL.lambda=lambda;
    
    PL.HRMS2=NaN;
    PL.VRMS2=NaN;
    PL.PRMS2=NaN;
    PL.HSlopeMax2=NaN;
    PL.VSlopeMax2=NaN;
    PL.HPL2=NaN;
    PL.VPL2=NaN;
    PL.lambda2=lambda;
else
    PRNindexforLocalTest=-1; % solution reliable - non � da escludere
    %% Protection Level 
    % Kai Borre computation
    % Conversion from ECEF to ENU coordinates
    xyzdtllh=ECEF2LLH(xyzdt(1:3));
    phi=xyzdtllh(1)*pi/180; %phi
    lambdAngle=xyzdtllh(2)*pi/180; %lambda
    cl = cos(lambdAngle);  sl = sin(lambdAngle);
    cb = cos(phi);  sb = sin(phi);
    F = [-sl -sb*cl cb*cl;
        cl -sb*sl cb*sl;
        0   cb      sb];
    H_ENU=[-RAIMparam.H(:,1:3)*F RAIMparam.H(:,4:end)];
    HWH_ENU=(H_ENU'*W*H_ENU)^(-1);
    
    %A=HWH*RAIMparam.H'*W;
    A=HWH_ENU*H_ENU'*W;
    
    PL.HRMS=sqrt(HWH_ENU(1,1)+HWH_ENU(2,2)); % Horizontal accuracy
    PL.VRMS=sqrt(HWH_ENU(3,3)); % Vertical accuracy
    PL.PRMS=sqrt(PL.HRMS^2+PL.VRMS^2);% Position accuracy
    % Calculated Bias vector
    B=H_ENU*A;    
    S=eye(size(B))-B;
    
    % Slopes vector
    HSlope=sqrt(diag(cov_noise)'.*(A(1,:).^2+A(2,:).^2)./diag(S)');
    VSlope=sqrt(diag(cov_noise)'.*A(3,:).^2./diag(S)');
    
    %HSlope=sqrt(A(1,:).^2+A(2,:).^2)./sqrt(diag(RAIMparam.S))';
    %VSlope=sqrt(A(3,:).^2)./sqrt(diag(S))';
    PL.HSlopeMax=max(HSlope);
    PL.VSlopeMax=max(VSlope);
    %trace(RAIMparam.S) % = Nsat-4
    PL.HPL=PL.HSlopeMax*sqrt(lambda);
    PL.VPL=PL.VSlopeMax*sqrt(lambda);
    PL.lambda=lambda;
    
    %%% Todd Walter computation
    T_thres=sqrt(chi2inv(1-settings.RAIM.PfaLocal,nu)); % sqrt(gammaincinv(1-settings.RAIM.PfaLocal,(nu-4)/2));
    k_Pmd=settings.RAIM.N_of_sd;
    
    K=HWH*RAIMparam.H'*W;
    PL.HRMS2=sqrt(HWH(1,1)+HWH(2,2)); % Horizontal accuracy
    PL.VRMS2=sqrt(HWH(3,3)); % Vertical accuracy
    PL.PRMS2=sqrt(PL.HRMS2^2+PL.VRMS2^2);% Position accuracy
    % Calculated Bias vector
    P=RAIMparam.H*K;     S=eye(size(P))-P;
    % Slopes vector
    HSlope=sqrt(diag(cov_noise)'.*(K(1,:).^2+K(2,:).^2)./diag(S)');
    VSlope=sqrt(diag(cov_noise)'.*K(3,:).^2./diag(S)');
    
    PL.HSlopeMax2=max(HSlope);
    PL.VSlopeMax2=max(VSlope);
    %trace(RAIMparam.S) % = Nsat-4
    PL.HPL2=PL.HSlopeMax2*T_thres+k_Pmd*PL.HRMS2;
    PL.VPL2=PL.VSlopeMax2*T_thres+k_Pmd*PL.VRMS2;
    PL.lambda2=T_thres;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If I want to avoid satellites exclusion(only detection and
% identification), enable following line:
%all_is_ok=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots...
% Plot distributions
% if 1==0
%     figure
%     xxbins=-8:0.01:8;
%     subplot(211)
%     plot(xxbins,normpdf(xxbins,0,1),xxbins,normpdf(xxbins,sqrt(lambda),1)), grid minor, hold on
%     plot(ones(size(0:.01:max(normpdf(xxbins,0,1))))*norminv(1-alpha0/2,0,1),0:.01:max(normpdf(xxbins,0,1)),'r-'), hold on, % alpha0/2 threshold
%
%     zeropos=ceil(length(xxbins)/2);
%     [minn, alpha0pos]=min(abs(xxbins-norminv(1-alpha0/2,0,1)));
%     alpha0check=1-2*trapz(xxbins(zeropos:alpha0pos),normpdf(xxbins(zeropos:alpha0pos),0,1))
%
%     subplot(212)
%     plot(xbins,ncx2pdf(xbins,nu,lambda),'r.-',xbins,chi2standard,'.-'), grid minor, hold on
%     plot(ones(size(0:.01:max(chi2standard)))*Globalthresh,0:.01:max(chi2standard),'k')
% end


end

