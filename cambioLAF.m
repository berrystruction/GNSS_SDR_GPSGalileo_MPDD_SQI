% %%
clc
close all,

timesize=length(trackResults(1).LAF);
steppy=1;
%settingsMPDD.mvgAvgtime=settingsMPDD.mvgAvgtime*100;
kkk=0;
wcoeff=zeros(settingsMPDD.M,timesize,length(trackResults));
U=zeros(length(trackResults(1).LAF(1).D)-kkk,settingsMPDD.M,timesize,length(trackResults));
D=zeros(length(trackResults(1).LAF(1).D)-kkk,timesize,length(trackResults));
y=D;

wcoeffk2=zeros(settingsMPDD.M,timesize);


trackResults2=trackResults;

channelNr=1;
h = waitbar(0);
for channelNr = 1:length(trackResults)
    
    x=0;
    set(h,'Name',['LAF: siamo al channel ' num2str(channelNr) ' di ' num2str(length(trackResults))]);
    
    for timeindex=1:steppy:timesize
        waitbar(x,h,[num2str(round(x*100)) '%'])
        x=timeindex/timesize;
        %wcoeff(:,timeindex,channelNr)=trackResults(channelNr).LAF(timeindex).w;
        %wcoeff(:,timeindex,channelNr)=trackResults(channelNr).LAF_CLS(timeindex).w;
        %wcoeff(:,timeindex,channelNr)=trackResults(channelNr).LAF(timeindex).w;
        %wcoeff(:,timeindex,channelNr)=trackResults(channelNr).LAF_CLS_N(timeindex).w;
        %wcoeff(:,timeindex,channelNr)=trackResults(channelNr).LAF(timeindex).w;
        %
        %     D(:,timeindex)=trackResults(channelNr).LAF(timeindex).D;
        %     U(:,:,timeindex)=trackResults(channelNr).LAF(timeindex).U;
        %     y(:,timeindex)=trackResults(channelNr).LAF(timeindex).y;
        
        if 1==0
            ausD=trackResults(channelNr).LAF(timeindex).D;
            d(:,timeindex)=ausD(ceil(length(ausD)/2)-settingsMPDD.multicorr.Npoints:ceil(length(ausD)/2)+settingsMPDD.multicorr.Npoints);
            
            ausU=trackResults(channelNr).LAF(timeindex).U;
            u(:,timeindex)=ausU(ceil(length(ausU)/2)-settingsMPDD.multicorr.Npoints:ceil(length(ausU)/2)+settingsMPDD.multicorr.Npoints,ceil(settingsMPDD.M/2));
            
            dreference=trackResults(channelNr).LAF(5).D;
            dreference=dreference(ceil(length(ausD)/2)-settingsMPDD.multicorr.Npoints:ceil(length(ausD)/2)+settingsMPDD.multicorr.Npoints);
        else
            %u(:,timeindex)=trackResults(channelNr).LAF(timeindex).U(9:33,9);
            u(:,timeindex)=trackResults(channelNr).LAF(15).d;
            %u(:,timeindex)=mean([trackResults(channelNr).LAF([1:15]).d],2);
            
            %d(:,timeindex)=trackResults(channelNr).LAF(timeindex).D(9:33);
            dreference=trackResults(channelNr).LAF(timeindex).d;%max(trackResults(channelNr).LAF(5).d);
            %dreference=mean([trackResults(channelNr).LAF(10:100).d],2);
        end
        
        Dmeasur=dreference;%d(:,timeindex);%/trapz(d(:,timeindex));%/max(d(:,timeindex))
        [LAF]=LAF_SQM(u(:,timeindex),Dmeasur,settingsMPDD.M,3);
        
        
        %D(:,timeindex,channelNr)=LAF.D;
        %y(:,timeindex,channelNr)=LAF.y;
        %U(:,1:settingsMPDD.M,timeindex,channelNr)=LAF.U;
        %         LAFaus.D=LAF.D;  LAFaus.N=LAF.N;  LAFaus.U=LAF.U;
        %         LAFaus.e=LAF.e;   LAFaus.y=LAF.y;  LAFaus.ynew=LAF.ynew; LAFaus.w=LAF.w;
        
        trackResults2(channelNr).LAF(timeindex)=LAF;
        trackResults2(channelNr).CN0=trackResults(channelNr).CN0;
    end
    
end




close(h);

%%
[SQIchannels, MPDDpoints]=MPDDcomputation(trackResults,settings,settingsMPDD,FLL_results);

[SQIchannels, MPDDpoints]=MPDDcomputation(trackResults2,settings,settingsMPDD,FLL_results);

[SQIchannels,MPDDpoints,LAFerrors,LAFerrorsBiased,Idealerror,ee]=SymmetricGaussianMethod(trackResults,settings,settingsMPDD,FLL_results);
%%
for channelNr = 1:length(trackResults)
    trackResults2(channelNr).CN0=trackResults(channelNr).CN0;
    trackResults2(channelNr).PRN=trackResults(channelNr).PRN; 
end
%%
[SQIchannels,MPDDpoints,LAFerrors,LAFerrorsBiased,Idealerror,ee]=SymmetricGaussianMethod(trackResults,settings,settingsMPDD,FLL_results);
for chn=1:length(trackResults), figure, plot(LAFerrors(chn,:),'bo-'), hold on, plot(LAFerrorsBiased(chn,:),'r.-'), 
    plot(Idealerror(chn,:),'k.-'),
    hold off,grid minor,title(['Channel ' num2str(chn)]), legend('Unbiased error','Biased error','Ideal error'),
end

% figure
% for kk=1:length(trackResults.LAF) 
%     plot(trackResults.LAF(kk).D,'o-'), grid on, title(num2str(settingsMPDD.mvgAvgtime*1e-3*kk)), pause(0.1);
% end

%%
% figure,
for channelNr = [1 5 6 7]
    for iii=1:length(trackResults(channelNr).LAF),
        %     subplot(311),
        %     plot(trackResults(channelNr).LAF(iii).D,'.-'),grid on, hold on, plot(trackResults2(channelNr).LAF(iii).u,'o-'),title(['Time: ' num2str(settings.seek_sec+settingsMPDD.mvgAvgtime*1e-3*iii)]),
        %     hold off,
        %     subplot(312),
        %     stem(trackResults(channelNr).LAF(iii).w),
        %     grid on, title('LS')
        %     subplot(313),    stem(trackResults2(channelNr).LAF(iii).w), grid on, title('CLS')
        %     pause(0.3),
        
        avgLSerror(channelNr,iii)=mean(trackResults(channelNr).LAF(iii).e(6:28));
        avgCLSerror(channelNr,iii)=mean(trackResults2(channelNr).LAF(iii).e(6:28));
        
    end
    
    figure, plot(avgLSerror(channelNr,:),'x-'), hold on, plot(avgCLSerror(channelNr,:),'o-'), grid on, title(['Approximation error of Channel ' num2str(channelNr)]), legend('LS','CLS')
end

figure,
iii=540;
plot(trackResults2(channelNr).LAF(iii).D,'o-'), hold on, grid on, plot(trackResults2(channelNr).LAF(iii).y,'*-')

