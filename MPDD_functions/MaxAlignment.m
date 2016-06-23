function [new_d]=MaxAlignment(d,newdim,indexMax) % indexMax
[~, dindexMax]=max(d);
difference=abs(dindexMax-indexMax);
%e=[zeros(difference,1); d ; zeros(difference,1)];
e=[d(1)*ones(difference,1); d ; d(end)*ones(difference,1);];


dindexMax=dindexMax+difference;
%dindexMax=ceil((length(d))/2);
%if length(d)> (dindexMax+(newdim-1)/2)
%cc=length(d)
% ee=[dindexMax-(newdim-1)/2, dindexMax ,dindexMax+(newdim-1)/2]

new_d=[e((dindexMax-(newdim-1)/2):1:(dindexMax+(newdim-1)/2))];
% for cc=1:M
%     new_d(cc)=new_d(M+1);
% end

% 
% 

end