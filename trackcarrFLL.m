function [FLLres]=trackcarrFLL(fid,acqResults,settings,PRN)

% FLLres.StatusFLL= zeros(1,32);
% FLLres.cnt_skp= zeros(1,32);
% FLLres.doppler_estFLL= zeros(1,32);
% FLLres.phaseFLL= zeros(1,32);
% FLLres.seeksec_corr= zeros(1,32);


%for PRN = settings.acqSatelliteList
    
    if acqResults.peakMetric(PRN) >= settings.acqThreshold %acqResults.carrFreq(PRN)~=0
        code_phase=acqResults.codePhase(PRN);
        doppler_est=acqResults.carrFreq(PRN);
        
        [StatusFLL,cnt_skp,doppler_estFLL,code_phaseFLL] = trackcarrFLL_DLL(fid,settings,code_phase,doppler_est,PRN);
        fprintf('IF Doppler [Hz]: %d \n',doppler_estFLL-settings.IF);
        
        %settings.seek_sec=settings.seek_sec+(code_phaseFLL/settings.samplingFreq-settings.seek_sec);
        FLLres.StatusFLL=StatusFLL;
        FLLres.cnt_skp=cnt_skp;
        FLLres.doppler_estFLL=doppler_estFLL;
        FLLres.phaseFLL=code_phaseFLL;
        FLLres.seeksec_corr=code_phaseFLL/settings.samplingFreq-settings.seek_sec;
    else
        FLLres.StatusFLL=0;
        FLLres.cnt_skp=-1;
        FLLres.doppler_estFLL=-1;
        FLLres.phaseFLL=-1;
        FLLres.seeksec_corr=-1;
        
    end
    
%end




end