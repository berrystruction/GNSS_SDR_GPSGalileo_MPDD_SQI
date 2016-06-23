function  [r, kmax]  = Hk(v,k)
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