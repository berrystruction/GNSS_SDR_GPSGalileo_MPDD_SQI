%function [ wnormalized,est_alpha,norms ] = w_normaliz(w,timesize)
% wnormalized=zeros(size(w));
% est_alpha=zeros(1,timesize);
% norms=est_alpha;
%
% for timeindex=1:timesize
%    wnormalized(:,timeindex)= w(:,timeindex)/max(w(:,timeindex));
%    est_alpha(timeindex)=sum(wnormalized(2:end,timeindex).^2);
%    norms(timeindex)=sum(wnormalized(1:end,timeindex).^2);
% end

function [ wk, wnorm ] = wHardTresh(w,timesize,k)

wk=zeros(size(w));
wnorm=wk;
for timeindex=1:timesize
    wnorm(:,timeindex)= (w(:,timeindex)/max(abs(w(:,timeindex))));
    %wk(:,timeindex)= HardThresh(abs(wnorm(:,timeindex)),k);
    wk(:,timeindex)=Hk((wnorm(:,timeindex)),k);
end


end

function  [r kmax]  = Hk(v,k)
% implementation of hard thresholding operator
vec=abs(v);
kmax=[];
minn=min(vec);
for indexkmax=1:k % i primi k elementi massimi in ampiezza
    [maxx indexmax]=max(vec);
    kmax=[kmax,indexmax];
    lenindexmax=length(indexmax);
    for Ntimesmax=1:lenindexmax % Ciclo se un massimo ha più occorrenze
        vec(indexmax(Ntimesmax))=minn-1;
    end
    
end
r=zeros(size(v));
r(kmax)=v(kmax);
end

function x = HardThresh(y,t)
% HardThresh -- Apply Hard Threshold 
%  Usage 
%    x = HardThresh(y,t)
%  Inputs 
%    y     Noisy Data 
%    t     Threshold
%  Outputs 
%    x     y 1_{|y|>t}
%
	x   = y .* (abs(y) > t);
end