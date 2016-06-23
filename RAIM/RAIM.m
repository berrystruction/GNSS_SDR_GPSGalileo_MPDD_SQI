function [RAIMoutput] = RAIM(dx,delta_rho,H,PAPrioriCovariance)
%RAIM 
%Convert the whole thing to ENU
% uposllh=ECEF2LLH(upos(1:3));
% phi=uposllh(1)*pi/180; %phi
% lambda=uposllh(2)*pi/180; %lambda
% cl = cos(lambda);  sl = sin(lambda);
% cb = cos(phi);  sb = sin(phi);
% F = [-sl -sb*cl cb*cl;
%     cl -sb*sl cb*sl;
%     0   cb      sb];
% H0=[-H(:,1:3)*F H(:,4:end)];
% A=(H0'*H0)\H0';


%parity vector, see if that works!
%[Q R]=qr(H);
%P=Q(:,5:end)'; % Q^T
%p=P*delta_rho;

%S=P'*P;

%S=eye(numsat)-H*(H'*H)^-1*H';
w=delta_rho-H*dx; % residual vector
%w=S*delta_rho;
%RAIMoutput.teststat=sqrt((w'*w)/(numsat-4));
%RAIMoutput.sig2=sqrt(sum(diag(w*w')));
%sig2=(w'*w)/(numsat-4);
% sig2Y=var(delta_rho);
% Cv=sig2*S;

% ppp=p'*p www=w'*w For test, these values are equals

%RAIMoutput.p=p;
%RAIMoutput.P=P;
RAIMoutput.w=w; RAIMoutput.y=delta_rho; RAIMoutput.yls=H*dx;
RAIMoutput.H=H;
RAIMoutput.aprioriPk=PAPrioriCovariance;
%RAIMoutput.S=S;
%RAIMoutput.SSE=w'*w;
%RAIMoutput.sig2=sig2;

%For the Protection Level
% A=(H'*H)^(-1)*H';
% HSlope=sqrt(A(1,:).^2+A(2,:).^2)./sqrt(diag(S))';
% VSlope=sqrt(A(3,:).^2)./sqrt(diag(S))';
% 
% HSlopeMax=max(HSlope);

% threshold=chi2inv(1-settings.RAIM.Pfa,numsat-4);

% navigation error
% est_e=w./diag(S);
% navErr=inv(H0'*H0)*H0'*est_e;
% w2=S*est_e;
% 
% pp=P*est_e;
% sum((p-pp).^2);

% %Now calculate the detection threshold
%thresraim=icdf('chi2',1-settings.RAIM.Pfa,numsat-4)/sqrt(numsat-4);
% % numsatall=numsat;
% 
% X=norminv(settings.RAIM.Pmd,0,1);
% l1=sqrt(20.5);l2=sqrt(30.5);
% 
% 
% %Protection level
% RAIMoutput.HSL=max(HSlope); RAIMoutput.VSL=max(VSlope);
% k1=chi2app(sqrt(thresraim),numsat-4,l1);
% k2=chi2app(sqrt(thresraim),numsat-4,l2);
% RAIMoutput.lambda_sqrt=l2-(k2-X)*(l2-l1)/(k2-k1);
% sig2p=sqrt(sum(diag(parity(:,:,step))));
% pbias=sig2p*lambda_sqrt;
% RAIMoutput.HPL=pbias*HSL;
% RAIMoutput.VPL=pbias*VSL;




end

