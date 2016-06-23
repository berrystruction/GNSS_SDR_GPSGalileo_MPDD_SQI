function [N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAFmethod_autocorrelation_IHT(u,d,M,LOS,Num_epoch,k)

N=length(u);
%M=floor(Particular_ratio*(N-1));
Particular_ratio=M/N;

rownumber=N+M-1; % N-M+1


U=zeros(rownumber,M);           % Correlation Matrix
w=-50*ones(M,Num_epoch);           % weights of ideal correlation
error_vec=zeros(M,Num_epoch);   % error vector for covariance calculation

powers=zeros(Num_epoch,5);  % Powers matrix

% LOS=1;                    % Index of LOS
% ynewMean=zeros(N-M+1,1);
% yMean=zeros(N-M+1,1);

y=zeros(rownumber,Num_epoch);
ynew=zeros(rownumber,Num_epoch);


e=zeros(rownumber,1);

% extension of correlations with the border points
for col=1:M
    U(1:floor(rownumber/2),col)=u(1)*ones(floor(rownumber/2),1);
    U(floor(rownumber/2)+2:end,col)=u(end)*ones(floor(rownumber/2),1);
    U(col:N+col-1,col)=u;%u(M-col+1:1:N-col+1);
end
%phi=U'*U;

epoch=1;

if length(d)~=N+M-1
    D=zeros(rownumber,1);
    D(1:floor(rownumber/2))=d(1)*ones(floor(rownumber/2),1);
    D(floor(rownumber/2)+2:end)=d(end)*ones(floor(rownumber/2),1);
    D(floor(rownumber/2)+1-floor(N/2):floor(rownumber/2)+1+floor(N/2))=d;
else
    D=d(:);
end
%theta=U'*D;
%w(:,epoch)=(phi)\theta; % Estimation of FIR coefficients


% U1=U';
% D1=D';
% 
% cvx_begin quiet
%     variable w(1,M);
%     minimize( norm(w,1) );
%     subject to
%     (w*U1-D1)>=0;
%     %ones(1,M)*w'==1.5
% cvx_end
% w=w';


Uorig=U;
Dorig=D;

Niteration=200;
U=U*0.97/norm(U);
D=D*0.97/norm(D);
for indexiter=1:Niteration
    w=Hk(w+U'*(D-U*w),k);
%     if(rem(indexiter,20))==0
%         w
%     end
end




%%%%%
%y=filter(w(:,epoch),1,u(M:N,epoch));
y(:,epoch)=U*w(:,epoch);
%e=d(:,epoch)-y;
e(:,epoch)=D-y(:,epoch);

% Power measurement
%ynew=(y-w(1,epoch)*u(:,epoch));
ynew(:,epoch)=y(:,epoch)-w(LOS,epoch)*U(:,LOS);
%ynew=U(:,2:end)*w(2:end,epoch); % Same thing as ynew

% ynewMean=ynewMean+ynew(:,epoch);
% yMean=yMean+y(:,epoch);

% powers= [Pu Pside Ptot Pmeasur Pe]
%                                               u(M:N,epoch)
powers(epoch,1)=sum(abs(w(LOS,epoch)*U(:,LOS)).^2);        % Pu
powers(epoch,2)=sum(abs(ynew(:,epoch)).^2);                % Pside
powers(epoch,3)=sum(abs(y(:,epoch)).^2);                   % Ptot

powers(epoch,4)=sum(abs(D).^2);                            % Pmeasur
powers(epoch,5)=sum(abs(e(:,epoch)).^2);                   % Pe

% P=U*inv(U'*U)*U'; % Projection operator (Haykin pag.498)
% emin=(eye(size(P))-P)*D;
% result=emin'*D

error_vec(:,epoch)=w(:,epoch)-mean(w(:,epoch));%-w(1,epoch);

% powers2(epoch)=mean(y(:,epoch).^2+(w(LOS,epoch)*U(:,LOS)).^2);

end
