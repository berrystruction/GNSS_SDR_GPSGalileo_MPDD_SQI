function [indices, minn, h7, h8, theorVectors_m, distancesV]=MPDistanceDetectorAlternate(alpha,M,wFinal2,timesize,timevect)
tic
%M=M-2;
distancesV=zeros(M,timesize);
totnumcode=2^M/2;
theorVectors_m=zeros(M,totnumcode);
ii=1;
numbertoConv=1;
indexx=1;
% alpha=0.5; % per regolare i livelli dei raggi di MP
while numbertoConv<2^M
    if(numbertoConv==2^ii)
        numbertoConv=numbertoConv+1;
        ii=ii+1;
    end
    aus=dec2bin(numbertoConv,M);
%     if aus(M)=='0'
%         numbertoConv=numbertoConv+1
%         aus=dec2bin(numbertoConv,M);
%     end       
    for num=1:M-1 % genero il sigolo codice
        theorVectors_m(M-num+1,indexx)=str2double(aus(num))*alpha;
    end 
    theorVectors_m(1,indexx)=str2double(aus(M));
    numbertoConv=numbertoConv+2;
    indexx=indexx+1;
end

%if M==8
    %wIdealCN0
    %theorVectors_m(:,1)=abs(wIdeal42(1:M));
    %wFinal2=abs(wFinal2);
%end

h7=figure;
endpoint=size(theorVectors_m,2);
%% calcolo distanze dal vettore teorico
for Mindex= 1: endpoint%2^M-M
    specialindex=M;
    while specialindex>0 && theorVectors_m(specialindex,Mindex)==0
        specialindex=specialindex-1;
    end
    
    for timeIndex=1:timesize
        % normalization
        normalizedW=wFinal2(1:M,timeIndex)/sqrt(sum(wFinal2(1:M,timeIndex).^2));%max(wFinal2(1:M,timeIndex))+1;
        %normalizedW=normalizedW/sqrt(sum(normalizedW.^2));
        distancesV(Mindex,timeIndex)= sqrt(sum((theorVectors_m(:,Mindex)-normalizedW).^2));
    end
    if Mindex~=1
        plot(timevect,distancesV(Mindex,:)), grid on, hold all %[colorr(Mindex) '.-']    
    else
        plot(timevect,distancesV(Mindex,:),'.-'), grid on, hold all %[colorr(Mindex) '.-']
    end
end
title(['Distances between theor. LAF vectors and Measured LAF vectors (\alpha=' num2str(alpha) ')']);
legend('LOS')%,'LOS+1st MP Ray','LOS+2nd MP Ray','1st+2nd','LOS+1st+2nd')

[minn, indices]=min(distancesV(1:endpoint,:),[],1);
toc
h8=figure;
plot(timevect,indices,'.')
grid on
xlabel('Time [s]')
title(['MP Detector results with \alpha=' num2str(alpha)])


end