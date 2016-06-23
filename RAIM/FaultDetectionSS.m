function [all_is_ok, PRNindexforLocalTest,TGlobalTest,TlocalTest,Globalthresh,PL]=FaultDetectionSS(pos,RAIMparam,settings,cov_noise,activeChnList_GPS_L1,satpos, pseudoranges1_GPS_L1, channels_GPS_L1, Doppler_GPS_L1,SQI)

% Fault detection and exclusion function

Nmeasurement=size(RAIMparam.w,1);
%nu=Nmeasurement-4;

W=(cov_noise)^(-1);
HWH=(RAIMparam.H'*W*RAIMparam.H)^(-1);

A0=HWH*RAIMparam.H'*W; % for Full set solution

% % % % % % % % % % % covw=cov_noise-RAIMparam.H*HWH*RAIMparam.H';
% % % % % % % % % % % if nargin==5
% % % % % % % % % % %     % Test parameters for the KF
% % % % % % % % % % %     d=RAIMparam.w;
% % % % % % % % % % %     S=RAIMparam.H*RAIMparam.aprioriPk(1:4,1:4)*RAIMparam.H'+cov_noise;
% % % % % % % % % % %     TGlobalTest=sqrt(d'*S^(-1)*d);
% % % % % % % % % % %
% % % % % % % % % % %     wnormalized=d./sqrt(diag(covw));
% % % % % % % % % % % else
% % % % % % % % % % %     % Test parameters for the LS
% % % % % % % % % % %     wnormalized=RAIMparam.w./sqrt(diag(covw));
% % % % % % % % % % %     TGlobalTest=RAIMparam.w'*W*RAIMparam.w;
% % % % % % % % % % % end
% % % % % % % % % % %
% % % % % % % % % % % [TlocalTest, PRNindexforLocalTest]=max(abs(wnormalized)); % observation with the largest abs value is tested
% % % % % % % % % % %

d=zeros(1,Nmeasurement);
for trials=1:Nmeasurement
    positions=zeros(1,Nmeasurement); positions(trials)=1; positions=(positions==0);    
    activeChnList_GPS_L1_aus=activeChnList_GPS_L1(positions);

    % [pos,satElev(activeChnList_GPS_L1),az,dop] =
    pos_i=leastSquarePos_Doppler(satpos(:,activeChnList_GPS_L1_aus), pseudoranges1_GPS_L1(activeChnList_GPS_L1_aus)', channels_GPS_L1(activeChnList_GPS_L1_aus), Doppler_GPS_L1(activeChnList_GPS_L1_aus), settings);
    
    % Solution separation discriminators
    d(trials)=norm(pos-pos_i,2);
    
    H_i=RAIMparam.H; H_i(trials,:)=zeros(1,4);
    A_i=(H_i'*W*H_i)^(-1)*H_i'*W;

    dP=(A_i-A0)*cov_noise*(A_i-A0)';
    [lambda_i,D_i]=(eig(dP(1:2,1:2))); % Horizontal component
end
% ... CONTINUARE

Pmd = settings.RAIM.Pmd;%0.005; %0.2
alpha0 = settings.RAIM.PfaLocal;%0.001; % 0.107
%alpha = 0.001; LA CALCOLO

% lambda=(norminv(1-alpha0/2,0,1)+norminv(1-Pmd,0,1))^2;
% 
% % Globalthresh calculation
% xbins=0:.1:50;
% chi2standard=chi2pdf(xbins,nu);
% diff=abs(ncx2pdf(xbins,nu,lambda)-chi2standard);
% [minn, minpos]=min(diff(2:end)); % escludo lo zero nell'origine
% Globalthresh=xbins(minpos);
% 
% all_is_ok=1; % Initialization
% 
% % Global test
% if TGlobalTest>Globalthresh
%     GlobalTestresults=1; % failed
%     % Performing local test...
%     if TlocalTest<=norminv(1-alpha0/2,0,1);
%         PRNindexforLocalTest=0; % solution unreliable
%         all_is_ok=0;
%     else
%         all_is_ok=-1;
%     end
%     PL.HRMS=NaN;
%     PL.VRMS=NaN;
%     PL.PRMS=NaN;
%     PL.HSlopeMax=NaN;
%     PL.VSlopeMax=NaN;
%     PL.HPL=NaN;
%     PL.VPL=NaN;
%     PL.lambda=lambda;
% else
%     PRNindexforLocalTest=-1; % solution reliable - non � da escludere
%     %% Protection Level
%     % Conversion from ECEF to ENU coordinates
%     xyzdtllh=ECEF2LLH(xyzdt(1:3));
%     phi=xyzdtllh(1)*pi/180; %phi
%     lambdAngle=xyzdtllh(2)*pi/180; %lambda
%     cl = cos(lambdAngle);  sl = sin(lambdAngle);
%     cb = cos(phi);  sb = sin(phi);
%     F = [-sl -sb*cl cb*cl;
%         cl -sb*sl cb*sl;
%         0   cb      sb];
%     H_ENU=[-RAIMparam.H(:,1:3)*F RAIMparam.H(:,4:end)];
%     HWH_ENU=(H_ENU'*W*H_ENU)^(-1);
%     
%     %A=HWH*RAIMparam.H'*W;
%     A=HWH_ENU*H_ENU'*W;
%     
%     PL.HRMS=sqrt(HWH_ENU(1,1)+HWH_ENU(2,2)); % Horizontal accuracy
%     PL.VRMS=sqrt(HWH_ENU(3,3)); % Vertical accuracy
%     PL.PRMS=sqrt(PL.HRMS^2+PL.VRMS^2);% Position accuracy
%     % Calculated Bias vector
%     B=H_ENU*A;
%     S=eye(size(B))-B;
%     
%     % Slopes vector
%     HSlope=sqrt(diag(cov_noise)'.*(A(1,:).^2+A(2,:).^2)./diag(S)');
%     VSlope=sqrt(diag(cov_noise)'.*A(3,:).^2./diag(S)');
%     
%     %HSlope=sqrt(A(1,:).^2+A(2,:).^2)./sqrt(diag(RAIMparam.S))';
%     %VSlope=sqrt(A(3,:).^2)./sqrt(diag(S))';
%     PL.HSlopeMax=max(HSlope);
%     PL.VSlopeMax=max(VSlope);
%     %trace(RAIMparam.S) % = Nsat-4
%     PL.HPL=PL.HSlopeMax*sqrt(lambda);
%     PL.VPL=PL.VSlopeMax*sqrt(lambda);
%     PL.lambda=lambda;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If I want to avoid satellites exclusion(only detection and
% identification), enable following line:
%all_is_ok=1;


end

