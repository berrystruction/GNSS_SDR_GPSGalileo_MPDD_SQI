%function [NewactiveChnList_GPS_L1,SVexcluded,SVexclIndex,unreliableSol,unrelindex]=ExclusionSatellites(currMeasNr,all_is_ok,navSolutions,activeChnList_GPS_L1,measurIndtoexclude,SVexclIndex,unrelindex)

% Exclusion
if all_is_ok<1
    
    if all_is_ok==-1 % SV exclusion from the active channel list
        
        SVexcluded(currMeasNr).PRN(SVexclIndex)=navSolutions.channel.PRN(activeChnList_GPS_L1(measurIndtoexclude),currMeasNr);
        SVexclIndex=SVexclIndex+1;
        
        if measurIndtoexclude==1 % se la misura � da escludere � la prima
            activeChnList_GPS_L1=activeChnList_GPS_L1(2:end);
            
        else
            
            if measurIndtoexclude==length(activeChnList_GPS_L1)
                activeChnList_GPS_L1=activeChnList_GPS_L1(1:end-1);
            else
                activeChnList_GPS_L1=activeChnList_GPS_L1([1:measurIndtoexclude-1 measurIndtoexclude+1:end]);
            end
        end
        
        if currMeasNr < startKalman
            pseudo=0; % da finire
            cov_noise=covNoiseMaker(channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,length(activeChnList_GPS_L1),pseudo,pseudo,SQI(activeChnList_GPS_L1,currMeasNr),SQIthreshold,settings.enableSQI);
        else
            KGPS = initSigmaKalmanGPS(KGPS,channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,kalmansettings.staticposition,pseudo,SQI(activeChnList_GPS_L1,currMeasNr),SQIthreshold,settings.enableSQI);
        end
    else
        
        if settings.enableSQI==true % SQI could find the fault and avoids solution unreliable
            
            [~, minpos]=min(SQI([channels_GPS_L1(activeChnList_GPS_L1).PRN],currMeasNr));
            additionalPRNexcluded=channels_GPS_L1(activeChnList_GPS_L1(minpos)).PRN; 
            
            SVexcluded(currMeasNr).PRN(SVexclIndex)=-2;
            SVexclIndex=SVexclIndex+1;
            SVexcluded(currMeasNr).PRN(SVexclIndex)=additionalPRNexcluded;
            SVexclIndex=SVexclIndex+1;
            
            activeChnList_GPS_L1=setdiff(activeChnList_GPS_L1,activeChnList_GPS_L1(minpos));
            all_is_ok=-1;
            %unrelindex=unrelindex-1;
            %unreliableSol(unrelindex)=[];
            
            if currMeasNr < startKalman
                pseudo=0; % da finire
                cov_noise=covNoiseMaker(channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,length(activeChnList_GPS_L1),pseudo,pseudo,SQI(activeChnList_GPS_L1,currMeasNr),SQIthreshold);
            else
                KGPS = initSigmaKalmanGPS(KGPS,channels_GPS_L1(activeChnList_GPS_L1),currMeasNr,kalmansettings.staticposition,pseudo,SQI(activeChnList_GPS_L1,currMeasNr),SQIthreshold);
            end
            
        else
            unreliableSol(unrelindex)=currMeasNr;
            unrelindex=unrelindex+1;
        end
        
        
        
    end
else
    SVexcluded(currMeasNr).PRN(SVexclIndex)=0;
end % if all_is_ok<1


%end
