function RAIM_Execution(TGlobalTest,TlocalTest,GlobalThres,unreliableSol,settings,navSolutions)
figure
subplot(211)
timeaus=([0:length(navSolutions.E(end,:))-1]+floor(settings.msToProcess/settings.navSolPeriod-length(navSolutions.E(end,:))))*settings.navSolPeriod/1e3;
%timeaus=[1:length(TGlobalTest)];

plot(timeaus,ones(size(1:length(TGlobalTest)))*GlobalThres,'--',timeaus,TGlobalTest,'.','MarkerSize',10,'LineWidth',2),
grid minor, title('Global Test results'),
legend('Global test threshold','Global test statistic')%,'Unreliable solution'),
xlabel('Time [s]')
% ,timeaus(unreliableSol),TGlobalTest(unreliableSol),'ko'

subplot(212)
plot(timeaus(TGlobalTest>GlobalThres),TlocalTest(TGlobalTest>GlobalThres),'.',...
     timeaus,ones(size(timeaus))*norminv(1-settings.RAIM.PfaLocal/2,0,1),'-.',...
     timeaus(unreliableSol),TlocalTest(unreliableSol),'ko','MarkerSize',10,'LineWidth',2)
grid minor,
legend('Local Test statistic','Local test threshold','Unreliable solution'),
xlabel('Time [s]')
title('Local Test results')
end
