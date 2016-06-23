function [ Meancount,Meancenter,gaussOrNot] = gaussianity(meanVec,nbins)
%% Gaussianity
%nbins=50;
alpha=0.05;
figure
[Meancount, Meancenter]=hist(meanVec,nbins);
Area=trapz(Meancenter,Meancount);
Meancount=Meancount/Area;

bar(Meancenter,Meancount)
grid on,
hold on, plot(Meancenter,normpdf(Meancenter,mean(meanVec),std(meanVec)),'r-')
title('Distribution Estimation')
legend('Histogram','Estimated Gaussian')


%% Test kologorov-smirnov x95
%testvariable = (meanVec-mean(meanVec))/(std(meanVec));%/sqrt(TotTrials));
testvariable=meanVec;
%gaussOrNot = kstest(testvariable);
gaussOrNot = chi2gof(testvariable,'Alpha',alpha,'NBins',nbins);%'Ctrs',Meancenter);


if gaussOrNot==1
    disp('This variable is not Gaussian')
else
    disp(['Gaussian parameter estimator at the ' num2str(alpha*100) '% significance level'])
end







end

