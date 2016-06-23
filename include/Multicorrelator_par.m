function [new_u, new_d]=Multicorrelator_par(localcode,tempdatacos,tempdatasin,codephasestep,remcodephase,blksize,minlag,stepp,maxlag)
ca=localcode;
k=ceil(maxlag+0.1);% +0.1 per approssimare sempre all'intero superiore anche nel caso di maxlag intero 1.0, 2.0 ecc.
ca=[ca(end-k+1:end) ca ca(1:1+k-1)];

%define index into prompt code vector
tcode=remcodephase:codephasestep:((blksize-1)*codephasestep+remcodephase);
tcode2=ceil(tcode)+k;
promptcode=ca(tcode2);
%%
% maxlags=20; % number of points
% maxlag=1.5;
% step=maxlag/(2*maxlags+1); % delay step
% minlag=step;
%%
p_p = sum(promptcode .* promptcode);
p_i = sum(promptcode .* tempdatasin);
p_q = sum(promptcode .* tempdatacos);
p_qi= sqrt(p_q.^2+p_i.^2);

% initialization
%indexlag=1;
lag=minlag:stepp:maxlag;
lenlag=length(lag);

e_i=zeros(1,lenlag);
e_q=e_i; e_p=e_q;

l_i=zeros(1,lenlag);
l_q=e_i; l_p=e_q;

pp_e=zeros(1,lenlag);
pp_l=pp_e;

pool=parpool ('local',4);      % Call to open the distributed processing
parfor indexlag=1:lenlag
    %define index into early code vector
    tcode=(remcodephase-lag(indexlag)):codephasestep:((blksize-1)*codephasestep+remcodephase-lag(indexlag));
    tcode2=ceil(tcode)+k;
    earlycode=ca(tcode2);
    
    %define index into late code vector
    tcode=(remcodephase+lag(indexlag)):codephasestep:((blksize-1)*codephasestep+remcodephase+lag(indexlag));
    tcode2=ceil(tcode)+k;
    latecode=ca(tcode2);
    
    %now get early, late, and prompt values for each
    e_i(indexlag) = sum(earlycode .* tempdatasin);
    e_q(indexlag) = sum(earlycode .* tempdatacos);
    e_p(indexlag) = sum(earlycode .* promptcode);
    
    l_i(indexlag) = sum(latecode .* tempdatasin);
    l_q(indexlag) = sum(latecode .* tempdatacos);
    l_p(indexlag) = sum(latecode .* promptcode);
    
    
    pp_e(indexlag)= sign(p_i)*sign(e_i(indexlag))*sqrt(e_i(indexlag).^2+e_q(indexlag).^2);
    pp_l(indexlag)= sign(p_i)*sign(l_i(indexlag))*sqrt(l_i(indexlag).^2+l_q(indexlag).^2);
    
    %indexlag=indexlag+1;
end

%% NO sign of the correlation
% new_u=[e_p(end:-1:1) p_p l_p]';
% pp_e= sqrt(e_i.^2+e_q.^2);
% pp_l= sqrt(l_i.^2+l_q.^2);
% new_d=[pp_e(end:-1:1) p_qi pp_l];

%% sign of the correlation saved
new_u=[e_p(end:-1:1) p_p l_p]';
new_d=[pp_e(end:-1:1) p_qi pp_l];


%tau=([-(length(new_d)-1)/2:(length(new_d)-1)/2]'*step)/codefreq;

new_d=new_d(end:-1:1)';
%new_u=new_u(end:-1:1);

delete(pool);           % Close the distributed computing
end