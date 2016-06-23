function [xCurrentEstimState,PAPostCovariance,KGain,PAPrioriCovariance] = KFroutine(zCurrentMeasurement,xPreviousEstimState,PAPostCovariance,HObservationMatrix,RObservationNoiseCovariance,PHIStateTransitMatrix,QStateNoiseCovariance,alpha)
% 
% Discrete Kalman filter - Core iteration
%
% ---
% Creation date: 11/09/2008
% Last mofication date: 
% --
% Author: EF
% State Prediction
%keyboard
xPredictedState = PHIStateTransitMatrix*xPreviousEstimState;

% A-priori error covariance matrix extrapolation
PAPrioriCovariance = alpha^2 * PHIStateTransitMatrix*PAPostCovariance*PHIStateTransitMatrix' + QStateNoiseCovariance;

% Kalman Gain Matrix
KGain = PAPrioriCovariance * HObservationMatrix' * inv(HObservationMatrix*PAPrioriCovariance*HObservationMatrix'+RObservationNoiseCovariance);

% State Estimate Update  
Kalmanresidual=zCurrentMeasurement - HObservationMatrix*xPredictedState;

xCurrentEstimState = xPredictedState + KGain*(Kalmanresidual);

% A-posteriori error covariance matrix
PAPostCovariance = (eye(length(xPreviousEstimState)) - KGain*HObservationMatrix) * PAPrioriCovariance * (eye(length(xPreviousEstimState)) - KGain*HObservationMatrix).' + KGain*RObservationNoiseCovariance*KGain.'; 
% to guarantee symmetry of the matrix  
end