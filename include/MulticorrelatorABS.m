function [new_u, new_d]=MulticorrelatorABS(localcode,tempdatacos,tempdatasin,codephasestep,remcodephase,blksize,minlag,stepp,maxlag)
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

vecindex=1;
for lag=minlag:stepp:maxlag
    %define index into early code vector
    tcode=(remcodephase-lag):codephasestep:((blksize-1)*codephasestep+remcodephase-lag);
    tcode2=ceil(tcode)+k;
    earlycode=ca(tcode2);
    
    %define index into late code vector
    tcode=(remcodephase+lag):codephasestep:((blksize-1)*codephasestep+remcodephase+lag);
    tcode2=ceil(tcode)+k;
    latecode=ca(tcode2);
    
    %now get early, late, and prompt values for each
    e_i(vecindex) = sum(earlycode .* tempdatasin);
    e_q(vecindex) = sum(earlycode .* tempdatacos);
    e_p(vecindex) = sum(earlycode .* promptcode);
    
    l_i(vecindex) = sum(latecode .* tempdatasin);
    l_q(vecindex) = sum(latecode .* tempdatacos);
    l_p(vecindex) = sum(latecode .* promptcode);
    
    vecindex=vecindex+1;
end

new_u=[e_p(end:-1:1) p_p l_p]';
pp_e= sqrt(e_i.^2+e_q.^2);
pp_l= sqrt(l_i.^2+l_q.^2);
new_d=[pp_e(end:-1:1) p_qi pp_l];

%tau=([-(length(new_d)-1)/2:(length(new_d)-1)/2]'*step)/codefreq;

new_d=new_d(end:-1:1)';
%new_u=new_u(end:-1:1);


end