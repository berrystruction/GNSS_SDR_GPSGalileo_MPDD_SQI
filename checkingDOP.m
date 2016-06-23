function   [newChnList,removedList]=checkingDOP(satPositions,rawP,satClkCorr,settings,TypeofSol,toRemoveList,activeChnList,NewactiveChnList)
% i.e. We fix DOP budget = 12
% other possibility is to create all possible combinations of channel
% list... not recommended
DOPbudget=settings.GDOPbudget;

[~,~,~,DOPactivechnlist,~]=leastSquarePos(satPositions(:,activeChnList), ...
    rawP(activeChnList) + satClkCorr(activeChnList) * settings.c, ...
    settings);

newChnList=NewactiveChnList;
%if TypeofSol==1
[~,~,~,DOP,~]=leastSquarePos(satPositions(:,newChnList), ...
    rawP(newChnList) + satClkCorr(newChnList) * settings.c, ...
    settings);
% else
%     %cov_noise=covNoiseMaker(trackResults(newChnList),currMeasNr,length(newChnList));
%     [~,~,~,DOP,WDOP,~] = WLS(satPositions(:,newChnList), ...
%         rawP(newChnList) + satClkCorr(newChnList) * settings.c, ...
%         settings,cov_noise);
% end

if DOP(1)<2*DOPactivechnlist(1)
    
    if DOP(1)>DOPbudget
        % Initialization
        PRNindex=1;
        DOP(1)=DOPbudget+1;
        newChnList=activeChnList;
        flagDOP=false;
        
        while (PRNindex<length(toRemoveList)) && flagDOP==false%(DOP(1)-DOPbudget)>0
            newChnList_prev=newChnList;
            newChnList=setdiff(activeChnList,toRemoveList(1:PRNindex));
            [~,~,~,DOP,~]=leastSquarePos(satPositions(:,newChnList), rawP(newChnList) + satClkCorr(newChnList) * settings.c,settings);
            if (DOP(1)>DOPbudget)
                flagDOP=true;
            end
            PRNindex=PRNindex+1;
        end
        
        if flagDOP==true
            newChnList=newChnList_prev;
        end
        
    end
    
    removedList=setdiff(activeChnList,newChnList);
    
    
else
    newChnList=activeChnList;
    removedList=[];
end


end