function SQI=SQIcomputation(trackResults,channelNr,MPresults,stima_CN0,settings,settingsMPDD,timesize,seekPRN)
%% SQI Part
%SQIchannels=zeros(32,floor(settings.msToProcess/settings.navSolPeriod));

if settingsMPDD.LAFversion==4 && exist('DistortionLAFindex','var')==1
    [ SQI ] = MPDD_part2(-MPresults,stima_CN0,settings,settingsMPDD,timesize,DistortionLAFindex);
else
    [ SQI ] = MPDD_part2(-MPresults,stima_CN0,settings,settingsMPDD,timesize);
end


%SQIchannels(trackResults(channelNr).PRN,:)=SQI;

% Print SQI results
figure
timevect=[settingsMPDD.mvgAvgtime:settingsMPDD.mvgAvgtime:timesize*settingsMPDD.mvgAvgtime]/1e3+(seekPRN); % Forse da correggere l asse dei tempi

subplot(311)
%     plot(timevect(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1),...
%         majorityIndex(settingsMPDD.MPDDoffset+floor(settings.navSolPeriod/settingsMPDD.mvgAvgtime)-1:floor(settings.msToProcess/settingsMPDD.mvgAvgtime)-1),'.','MarkerSize',12),
plot(timevect,MPresults,'.','MarkerSize',12),
grid minor
title(['MPDD output for Channel ' num2str(channelNr) ', PRN ' num2str(trackResults(channelNr).PRN)],'interpreter','Latex')
%     title(['Majority Rule for Channel ' num2str(channelNr) ', PRN ' num2str(trackResults(channelNr).PRN)],'interpreter','Latex')
legend('1: MP / -1: NO MP')

subplot(312)
austime=(settingsMPDD.MPDDoffset+[settings.navSolPeriod:settings.navSolPeriod:settings.msToProcess])/1e3+seekPRN; % asse dei tempi
plot(austime,SQI,'ro')
%ylim([-1.5 1.5])
grid minor
title(['SQI for PRN ' num2str(trackResults(channelNr).PRN) ' (with initial offset cutted)'],'interpreter','Latex')
xlabel('Time [s]')
ylim([0 1]);

subplot(313)
plot(timevect,stima_CN0(1:length(timevect)),'.-')
grid on
title('C/N_0')
xlabel('Time [s]')

end