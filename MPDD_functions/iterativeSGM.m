function [SQI,Idealerror,MPresults,delayestimation,residualApprox]=iterativeSGM(trackResults,settings,settingsMPDD,FLL_results,alpha)
% Gaussian parameters
% alpha = Significance level
% Pfa
% Pfa=10^-5;
% GaussThreshold=norminv(Pfa,0,sigma);

timesize=length(trackResults(1).LAF);
SQI=zeros(32,floor(settings.msToProcess/settings.navSolPeriod));

%LAFerrors
Idealerror=zeros(settings.numberOfChannels,timesize);
energyEvaluator=Idealerror;
availableComponents=energyEvaluator;
usedComponents=availableComponents;
MPresults=usedComponents;
pvalue=MPresults;

CalibrationTime=7; % [s]
CalibrationTimeIndex=(CalibrationTime/settingsMPDD.mvgAvgtime/1e-3);


lendiff=0;%length(trackResults(1).LAF(1).D)-length(trackResults(1).LAF(1).d);

NpointsLR=floor(length(trackResults(1).LAF(1).D(lendiff/2+1:end-lendiff/2))/2);

%Startlen=length(trackResults(1).LAF(1).e(1:end)); NpointsLR=(Startlen-1)/2;

first_coeff=1; end_coeff=length(trackResults(1).LAF(1).w)-0;
windexListM=zeros(length(trackResults(1).LAF(1).w(first_coeff:end_coeff))-1,timesize);

residualApprox=zeros(settings.numberOfChannels,2*NpointsLR+1,timesize);
delayestimation=zeros(settings.numberOfChannels,timesize);

for channelNr = 1:settings.numberOfChannels
    tic
    stima_CN0=trackResults(channelNr).CN0.CNo_SNV;
    %     timesize=length(trackResults(channelNr).LAF);
    %indice=1;
    %cc=zeros(len,1);
    err_vector=zeros(2*NpointsLR,timesize);
    
    for timeindex=1:timesize
        
        w_temp=trackResults(channelNr).LAF(timeindex).w(first_coeff:end_coeff);
        [~,windexList]=sort(abs(w_temp),'descend');
        % Choice of the w components - energy evaluation 90%
        corrEnergy=sum(trackResults(channelNr).LAF(timeindex).y.^2);
        lenw_temp=length(w_temp);
        
        %         kk_w=1;
        %         energy_temp=0;
        %         while  kk_w<=lenw_temp && energyEvaluator(channelNr,timeindex)<EnergyContent
        %             %windexList(kk_w)
        %             energy_temp=energy_temp+sign(w_temp(windexList(kk_w)))*sum((trackResults(channelNr).LAF(timeindex).U(:,windexList(kk_w))*w_temp(windexList(kk_w))).^2);%*sign(w_temp(kk_w));
        %             kk_w=kk_w+1;
        %             energyEvaluator(channelNr,timeindex)=energy_temp/corrEnergy*100;
        %         end
        
        windexList=setdiff(windexList,ceil(lenw_temp/2),'stable'); % remove LOS tap
        windexListM(:,timeindex)=windexList;
        
        %availableComponents(channelNr,timeindex)=kk_w-1-1;
        % Linea fake di controllo
        availableComponents(channelNr,timeindex)=length(windexList)+1;
        
        % Chi-square Testing part
        kk2_w=0;
        h=1;
        ccU=trackResults(channelNr).LAF(timeindex).U(:,first_coeff:end_coeff);
        while h==1 && kk2_w<availableComponents(channelNr,timeindex)
            kk2_w=kk2_w+1;
            
            if kk2_w>1 % first cicle: removing of LOS component
                approximErr=trackResults(channelNr).LAF(timeindex).D(lendiff/2+1:end-lendiff/2)...
                    -ccU(lendiff/2+1:end-lendiff/2,windexList(1:kk2_w-1))*w_temp(windexList(1:kk2_w-1));
                %approximErr=trackResults(channelNr).LAF(timeindex).D-trackResults(channelNr).LAF(timeindex).U(:,windexList(kk2_w))*w_temp(windexList(kk2_w));
            else
                approximErr=trackResults(channelNr).LAF(timeindex).D(lendiff/2+1:end-lendiff/2)...
                    -ccU(lendiff/2+1:end-lendiff/2,ceil(lenw_temp/2))*w_temp(ceil(lenw_temp/2));
                
                residualApprox(channelNr,:,timeindex)=approximErr;
                %approximErr=trackResults(channelNr).LAF(timeindex).D-trackResults(channelNr).LAF(timeindex).U(:,ceil(lenw_temp/2))*w_temp(ceil(lenw_temp/2));
            end
            %            Idealerror(channelNr,timeindex)=mean([approximErr([1:floor(len/2)]).^2; -(approximErr([floor(len/2)+2:end]).^2)]);
            err_vector(:,timeindex)=[approximErr(1:NpointsLR); -approximErr(NpointsLR+2:end)]';
            Idealerror(channelNr,timeindex)=mean(err_vector(:,timeindex));
            
            %plot(approximErr,'o-'),
            %hold on, plot([approximErr(1:NpointsLR); -approximErr(NpointsLR+2:end)],'*-'),grid on,
            %pause(.2)%, hold off
            if timeindex<CalibrationTimeIndex
                %[h] = chi2gof(Idealerror(channelNr,timeindex),'Alpha',alpha,'nbins',50);
                %[h] = kstest(Idealerror(channelNr,timeindex),'Alpha',alpha);%,'nbins',50);
                h=NaN;
                endCalibrationIndex=timeindex;
                sigma_err=std(Idealerror(channelNr,1:endCalibrationIndex),1);
                avg_err=mean(Idealerror(channelNr,1:endCalibrationIndex));
            else
                Idealerror(channelNr,timeindex)=(Idealerror(channelNr,timeindex)-avg_err)/sigma_err;
                %[h] = chi2gof(Idealerror(channelNr,timeindex),'Alpha',alpha,'nbins',50);
                [h,pvalue(channelNr,timeindex)] = kstest(Idealerror(channelNr,timeindex),'Alpha',alpha);
                
            end
        end
        
        usedComponents(channelNr,timeindex)=kk2_w;
        MPresults(channelNr,timeindex)=h;
        
        
        if usedComponents(channelNr,timeindex)>1 && h==1
            [~,posmax]=max(abs(residualApprox(channelNr,:,timeindex)));
            peakpos=ceil(length(residualApprox(channelNr,:,timeindex))/2);
            if posmax==peakpos
                [~,posmaxList]=sort(abs(residualApprox(channelNr,:,timeindex)),'descend');
                posmax=posmaxList(2);
            end
            delayestimation(channelNr,timeindex)=settingsMPDD.multicorr.multicorr_resolution*abs(posmax-peakpos);
            %gn=delayestimation(channelNr,timeindex); gn
        else
            delayestimation(channelNr,timeindex)=0;
        end
        
        kk_w=0;
        energy_temp=0;
        while  kk_w<=kk2_w-1
            kk_w=kk_w+1;
            %ccU=trackResults(channelNr).LAF(timeindex).U(:,first_coeff:end_coeff);
            if kk_w>1 % first cicle: removing of LOS component
                energy_temp=energy_temp+sign(w_temp(windexList(kk_w-1)))*sum((ccU(lendiff/2+1:end-lendiff/2,windexList(kk_w-1))*w_temp(windexList(kk_w-1))).^2);
            else
                energy_temp=energy_temp+sign(w_temp(ceil(lenw_temp/2)))*sum((ccU(lendiff/2+1:end-lendiff/2,ceil(lenw_temp/2))*w_temp(ceil(lenw_temp/2))).^2);
            end
            energyEvaluator(channelNr,timeindex)=energy_temp/corrEnergy*100;
        end
        
    end
    
    
    
    seekPRN=settings.seek_sec+FLL_results.seeksec_corr(trackResults(channelNr).PRN);%+ round((20-FLL_results.cnt_skp((trackResults(channelNr).PRN)))*settings.samplingFreq*1e-3);
    settingsMPDD.MPDDoffset=ceil((((ceil(settings.seek_sec)-settings.seek_sec)/settings.navSolPeriod*1e3)*1e3)/settingsMPDD.mvgAvgtime)+1;
    %settingsMPDD.MPDDoffset=ceil((((ceil(seekPRN)...
    %                        -(seekPRN))/settings.navSolPeriod*1e3)*1e3)/settingsMPDD.mvgAvgtime)+1;
    
    timevect=[settingsMPDD.mvgAvgtime:settingsMPDD.mvgAvgtime:timesize*settingsMPDD.mvgAvgtime]/1e3+(seekPRN);
    
    %MPresults(channelNr,:)=(MPresults(channelNr,:)*2-1);
    Idealerror(channelNr,1:CalibrationTimeIndex)=NaN;
    
    figure,
    subplot(321)
    plot(timevect,MPresults(channelNr,:),'bo'), grid minor, xlabel('Time [s]'), title(['Detection results for PRN ' num2str(trackResults(channelNr).PRN)]), legend('1 MP / 0 No MP')
    subplot(322),
    plot(timevect,Idealerror(channelNr,:),'ro'), grid minor, xlabel('Time [s]'), title('mean Residual Correlation Error')
    subplot(323)
    plot(timevect,delayestimation(channelNr,:)*settings.c), grid minor, xlabel('Time [s]'), ylabel('MP Delay [m]'), title('Delay estimation') %title('Energy approximated')
    subplot(324)
    plot(timevect,availableComponents(channelNr,:),'b',timevect,usedComponents(channelNr,:),'r'), grid minor, xlabel('Time [s]'), title('available vs used components'), legend('available','used')
    subplot(325)
    plot(timevect,stima_CN0(1:length(timevect)),'.-'), grid on, title('C/N_0 trend'), xlabel('Time [s]'), ylabel('C/N_0 [dBHz]')
    
    SQI(trackResults(channelNr).PRN,:)= MPDD_part2(-(MPresults(channelNr,:)*2-1),stima_CN0,settings,settingsMPDD,timesize);
    austime=(settingsMPDD.MPDDoffset+[settings.navSolPeriod:settings.navSolPeriod:settings.msToProcess])/1e3+seekPRN; % asse dei tempi
    ausindex=((CalibrationTime+austime(1))-austime)>1;    
    SQI(trackResults(channelNr).PRN,ausindex)=1; % NaN

    subplot(326)
    plot(austime, SQI(trackResults(channelNr).PRN,:),'ro'), grid on, title(['SQI for PRN ' num2str(trackResults(channelNr).PRN)])%,'interpreter','Latex')
    xlabel('Time [s]')
    ylim([0 1]);

    toc
    
    
    % Some plots
     if  channelNr==9 && 1==0
         figure, chn=channelNr; 
        for kk=1:5:timesize-40
            
            ccU=trackResults(chn).LAF(kk).U(:,first_coeff:end_coeff);
            if usedComponents(chn,kk)>1
                ccU=ccU(lendiff/2+1:end-lendiff/2,windexListM(1:usedComponents(chn,kk)-1,kk));
                
                usedComponents(chn,kk)
                ccw_aus=trackResults(chn).LAF(kk).w(first_coeff:end_coeff);
                ccw=ccw_aus(windexListM(1:usedComponents(chn,kk)-1,kk));
                ccUU=trackResults(chn).LAF(kk).U(lendiff/2+1:end-lendiff/2,ceil(lenw_temp/2))*ccw_aus(ceil(lenw_temp/2));
            else                
                ccw_aus=trackResults(chn).LAF(kk).w(first_coeff:end_coeff);
                ccw=zeros(length(trackResults(chn).LAF(kk).w(first_coeff:end_coeff))-1,1);
                ccUU=trackResults(chn).LAF(kk).U(lendiff/2+1:end-lendiff/2,ceil(lenw_temp/2))*ccw_aus(ceil(lenw_temp/2));
                ccU=trackResults(chn).LAF(kk).U(:,first_coeff:end_coeff-1);
            end
            subplot(221),
            plot(trackResults(chn).LAF(kk).D(lendiff/2+1:end-lendiff/2),'o-'), hold on,...
                plot(ccU*ccw+ccUU,'b+-'), hold on, plot(trackResults(chn).LAF(kk).y(lendiff/2+1:end-lendiff/2),'*-'), hold off
            title(['Time ' num2str(kk*settingsMPDD.mvgAvgtime*1e-3+settings.seek_sec) ' PRN ' num2str(trackResults(chn).PRN) ' timeindex ' num2str(kk) ]), grid on, legend('D','approx.','y')
            
            subplot(222),
            plot(ccU*ccw+ccUU,'b+-'), grid on, hold on, plot(ccUU,'*-'), grid on, hold on, plot(ccU*ccw,'.-'), grid on, legend('LOS+NLOS comp.','LOS comp.','NLOS comp.'), hold off
            subplot(223), stem(trackResults(chn).LAF(kk).w(first_coeff:end_coeff)), grid on, title(['# of used comp. ' num2str(usedComponents(chn,kk))])
            
            recomputederrdiff=(trackResults(chn).LAF(kk).D(lendiff/2+1:end-lendiff/2)-(ccU*ccw+ccUU));
            recomputeddiff=[recomputederrdiff(1:NpointsLR); -recomputederrdiff(NpointsLR+2:end)];
            
            subplot(224), plot(err_vector(:,kk),'*-'), grid on, hold on, plot(recomputeddiff,'^-'), hold off, ylim([-3e4 4e4])
            %plot(trackResults(chn).LAF(kk).D(lendiff/2+1:end-lendiff/2+1)-trackResults(chn).LAF(kk).y(lendiff/2+1:end-lendiff/2+1),'o-'), grid on
            pause(.4)
        end
    end
    
    
end




