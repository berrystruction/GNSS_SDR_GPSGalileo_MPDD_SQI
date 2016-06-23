function acqResults=myacquisition(data,settings)

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

fprintf('(');
for PRN = settings.acqSatelliteList
    
    
    [doppler_est, code_phase, status, ~,peakMetric] = signal_acquisition_FFT_RAPID(data,settings,PRN);
    
    if (status == 0) % If not acquired, retry with the complete Doppler range
        %         [doppler_est, code_phase, status, ~,peakMetric] = signal_acquisition_FFT(fid,T_int,settings,PRN);
        %
        %         if (status == 0)
        %             disp('no satellite...')
        %             return
        %         else
        acqResults.carrFreq(PRN)  = doppler_est;
        acqResults.codePhase(PRN) = code_phase;
        acqResults.peakMetric(PRN)= peakMetric;
        %         end
        fprintf('. ');

        
    else
        acqResults.carrFreq(PRN)  = doppler_est;
        acqResults.codePhase(PRN) = code_phase;
        acqResults.peakMetric(PRN)= peakMetric;
        fprintf('%02d ', PRN);

    end
    
    
end

fprintf(')\n');

acqResults.acqThreshold=settings.acqThreshold;
end

