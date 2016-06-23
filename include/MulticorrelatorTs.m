function [new_u, new_d, tau]=MulticorrelatorTs(promptcode,tempdatacos,tempdatasin,maxlags,offsetXcorrpoints,step,remcodephase,codefreq)

N=2*maxlags+1;
new_u=zeros(N,1);
new_d=new_u;
pp_i=new_d;
pp_q=pp_i;

code_m =codematrixgen(promptcode,maxlags,offsetXcorrpoints);


for ii = 1:N
    new_u(ii) = promptcode*code_m(:,ii);
end

for ii = 1:N
    pp_i(ii) = tempdatasin*code_m(:,ii);
end

for ii = 1:N
    pp_q(ii) = tempdatacos*code_m(:,ii);
end

new_d=sqrt(pp_i.^2+pp_q.^2);
tau=([-(length(new_d)-1)/2:(length(new_d)-1)/2]'*step+remcodephase)/codefreq;

end