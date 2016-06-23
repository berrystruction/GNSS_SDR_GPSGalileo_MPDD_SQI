%function [indices, minn, h7, theorVectors_m, distancesV]=MPDistanceDetector(alpha,M,wFinal2,timesize,timevect)
function [indices, minn, h7, theorVectors_m, distancesV, detectprofile, judmentfunc]=MPDistanceDetector(alpha,M,wFinal2,timesize,timevect,CN0profile,settings,settingsMPDD,PRN)

detectprofile=[]; 
judmentfunc=[];

%tic
%M=M-2;

distancesV=zeros(M,timesize);
totnumcode=2^floor(M/2);
theorVectors_m=zeros(M,totnumcode);
ii=1;
numbertoConv=1;
indexx=1;
% nuova generazione
Mnew=M;
M=ceil(M/2);

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
        theorVectors_m(ceil(Mnew/2)+num,indexx)= str2double(aus(num))*alpha;
    end
    theorVectors_m(ceil(Mnew/2),indexx)=str2double(aus(M));
    numbertoConv=numbertoConv+2;
    indexx=indexx+1;
end


mixedterms=theorVectors_m(:,1:totnumcode);
mixedterms(1:Mnew-M,:)=mixedterms(Mnew:-1:Mnew-M+2,:);
theorVectors_m=[theorVectors_m theorVectors_m(end:-1:1,:) mixedterms];

M=Mnew;


if settingsMPDD.plotMPDD==1 % if Plot is enabled
    h7=figure;
    subplot(221)
else
    h7='-';
end
endpoint=size(theorVectors_m,2);
%% calcolo distanze dal vettore teorico
for Mindex= 1: endpoint%2^M-M
    for timeIndex=1:timesize
        w=wFinal2(1:M,timeIndex);%/max(wFinal2(1:M,timeIndex));
        %distancesV(Mindex,timeIndex)=sqrt(sum((theorVectors_m(:,Mindex)-w).^2)); % distanza
        distancesV(Mindex,timeIndex)= sum((theorVectors_m(:,Mindex)-w).^2); % distanza^2
        %distancesV(Mindex,timeIndex)=sum(abs(theorVectors_m(:,Mindex)-w)); % norm L1 - Manhattan distance
    end
    
    if settingsMPDD.plotMPDD==1 % if Plot is enabled
        if Mindex~=1
            plot(timevect,distancesV(Mindex,:),'-'), grid on, hold on %[colorr(Mindex) '.-']
        else
            plot(timevect,distancesV(Mindex,:),'p-'), grid on, hold on %[colorr(Mindex) '.-']
        end
    end
end

if settingsMPDD.plotMPDD==1
    title(['Distances between theor. LAF vectors and Measured LAF vectors ($\alpha=' num2str(alpha) '$)'],'interpreter','Latex');
    legend('LOS');
    xlabel('Time [s]')
end

[minn, indices]=min(distancesV(1:endpoint,:),[],1);



%h8=figure;
if settingsMPDD.plotMPDD==1
    [ detectprofile, judmentfunc ] = MPDD_part2(2*(indices>1)-1,CN0profile,settings,settingsMPDD);
    
    subplot(223)
    austime=(settingsMPDD.MPDDoffset+[settings.navSolPeriod:settings.navSolPeriod:settings.msToProcess])/1e3+settings.seek_sec; % asse dei tempi
    %austime=timevect(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime):floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1);
    %austime=austime(1:end-1);%size(austime)
    
    plot(timevect,indices,'o')
    grid on
    xlabel('Time [s]')
    title(['MP Detector results with $\alpha=' num2str(alpha) '$' ],'interpreter','Latex')
    
    subplot(222)
    indicesNorm=2*(indices>1)-1;
    plot(timevect(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1),...
        indicesNorm(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1),'.',...
        austime,detectprofile,'o')
    %ylim([-1.5 1.5])
    grid minor
    title(['MPDD advanced (with initial offset cutted) for PRN ' num2str(PRN)],'interpreter','Latex')
    legend('MPDD','MPDD advanced')
    xlabel('Time [s]')
    
    subplot(224)
    plot(timevect,CN0profile)
    grid on
    title('C/N0','interpreter','Latex')
    xlabel('Time [s]')
    ylabel('dBHz')
%else
%    indices=2*(indices>1)-1;
end

%toc

end