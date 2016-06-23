function  [kmin,worstIndexList]  = SQImin(v,k,threshold)
% implementation of hard thresholding operator
vec=v;
kmin=[];
%minn=min(vec);
for indexkmin=1:k % i primi k elementi minimi in ampiezza
    [minn, indexmin]=min(vec);
    if minn<threshold
        kmin=[kmin,indexmin];
        vec(indexmin)=2;
    end
    %    lenindexmin=length(indexmin);
    %     for Ntimesmax=1:lenindexmin % Ciclo se un min ha più occorrenze
    %         vec(indexmin(Ntimesmax))=minn-1;
    %     end
    
end
worstIndexList=kmin;
%r=zeros(size(v));
%r(kmin)=v(kmin);
kmin=vec>1;
end