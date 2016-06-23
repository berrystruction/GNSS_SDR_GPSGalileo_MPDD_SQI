function [N,U,w,y,ynew,D,e,powers,error_vec,Particular_ratio]=LAF_New_method(u,d,M,LOS,tott)

for col=1:floor(tott/2)
    u=[u(1);u;u(end)];
    d=[d(1);d;d(end)];
end
% u=u/max(u);
% d=u/max(d);

% u(1:M)=u(M); u(end-M+1:end)=u(end-M+1);
% d(1:M)=d(M); d(end-M+1:end)=d(end-M+1);


N=length(u);
%M=floor(Particular_ratio*(N-1));
Particular_ratio=M/N;

U=zeros(N-M+1,M);           % Correlation Matrix
w=zeros(M,1);         % weights of ideal correlation
error_vec=zeros(M,1);   % error vector for covariance calculation

powers=zeros(1,5);  % Powers matrix

y=zeros(N-M+1,1);
ynew=zeros(N-M+1,1);


e=zeros(N-M+1,1);


for col=1:M
    U(:,col)=u(M-col+1:1:N-col+1);
    %U(:,M-col+1)=u(M-col+1:1:N-col+1);
end
%phi=U'*U;


if length(d)~=N-M+1
    D=d(M:N);
    %D=d(N:-1:M,epoch);
else
    D=d(:);
end
%theta=U'*D;
%w(:,1)=(phi)\theta; % Estimation of FIR coefficients
U1=U';
D1=D';

cvx_begin quiet
    variable w(1,M);
    minimize( norm(w*U1-D1) );
    subject to  
    %ones(1,M)*w'==1.5    
    w>=0;
cvx_end


w=w';

%%%%%
%y=filter(w(:,1),1,u(M:N,1));
y(:,1)=U*w(:,1);
%e=d(:,1)-y;
e(:,1)=D-y(:,1);

% Power measurement
%ynew=(y-w(1,1)*u(:,1));
ynew(:,1)=y(:,1)-w(LOS,1)*U(:,LOS);
%ynew=U(:,2:end)*w(2:end,1); % Same thing as ynew

% ynewMean=ynewMean+ynew(:,1);
% yMean=yMean+y(:,1);

% powers= [Pu Pside Ptot Pmeasur Pe]
%                                               u(M:N,1)
powers(1,1)=sum(abs(w(LOS,1)*U(:,LOS)).^2);        % Pu
powers(1,2)=sum(abs(ynew(:,1)).^2);                % Pside
powers(1,3)=sum(abs(y(:,1)).^2);                   % Ptot

powers(1,4)=sum(abs(D).^2);                            % Pmeasur
powers(1,5)=sum(abs(e(:,1)).^2);                   % Pe

% P=U*inv(U'*U)*U'; % Projection operator (Haykin pag.498)
% emin=(eye(size(P))-P)*D;
% result=emin'*D

error_vec(:,1)=w(:,1)-mean(w(:,1));%-w(1,1);





end