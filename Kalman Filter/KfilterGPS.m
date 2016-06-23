function [stateAPosteriori,PAPostCovariance,dop,RAIMresult] = KfilterGPS(K, pdr, dop, pos, vel, satPos, satVel, dt, ddt, AAA)

% PURPOSE:
% Computes vectors needed for Kalman filtering.
% -------------------------------------------------------------------------
% INPUTS:
% K             Kalman filter structure
% pdr           1st column satId 2nd pdr 3rd doppler
% pos           predicted position
% vel           predicted velocity
% satPositions  rows: x y z. columns: satId
% satClkCorr    clock corrections
% satVelocity   rows: x y z. columns: satId
% -------------------------------------------------------------------------
% OUTPUTS:
% stateAPosteriori
% PAPostCovariance
% -------------------------------------------------------------------------
% Created by: M. Rao
% Last update: M. Rao
% Version 1.0
% 27/07/2009
% -------------------------------------------------------------------------

noSat = numel(pdr);

Range = zeros(noSat,1);
DeltaRange = zeros(noSat,1);
U = zeros(noSat,4);
%keyboard
for IndexU = 1 : noSat
    %% Compute the Range with the relativistic correction and the satellite clock bias
    dX 		= (satPos(1,IndexU)*cos(AAA(IndexU)) + satPos(2,IndexU)*sin(AAA(IndexU))) - pos(1);
    dY 		= (satPos(2,IndexU)*cos(AAA(IndexU)) - satPos(1,IndexU)*sin(AAA(IndexU))) - pos(2);
    dZ 		= satPos(3,IndexU) - pos(3);
    
    dXv     = satVel(1,IndexU) - vel(1);
    dYv     = satVel(2,IndexU) - vel(2);
    dZv     = satVel(3,IndexU) - vel(3);
    
    Range(IndexU) = sqrt(power(dX,2)+power(dY,2)+power(dZ,2));
    DeltaRange(IndexU) = ((dXv*dX)+(dYv*dY)+(dZv*dZ))/Range(IndexU);
    U(IndexU,1:3) = [dX dY dZ]/norm(Range(IndexU));
    
%     HgeometricMatrix(IndexU,1)=satPos(1,IndexU)*cos(AAA(IndexU)) + satPos(2,IndexU)*sin(AAA(IndexU));
%     HgeometricMatrix(IndexU,2)=satPos(2,IndexU)*cos(AAA(IndexU)) - satPos(1,IndexU)*sin(AAA(IndexU));
%     HgeometricMatrix(IndexU,3)=satPos(3,IndexU);
%     HgeometricMatrix(IndexU,4)=1;
    
end
U(:,end) = -ones(noSat,1);
PrCorr = pdr-dt;
DeltaPrRaw = dop;
DeltaPrCorr = -DeltaPrRaw + ddt;

zCurrentMeasurement = [PrCorr; DeltaPrCorr];
NominalMeasurement = [Range; DeltaRange];
KFMeasurement = zCurrentMeasurement - NominalMeasurement;
HObservationMatrix = kron([1 0; 0 1],-U);
% stateAPosteriori = HObservationMatrix \ KFMeasurement;
% PAPostCovariance = K.PAPostCovariance;
[stateAPosteriori,PAPostCovariance,~,PAPrioriCovariance] = KFroutine(KFMeasurement,K.stateAPosteriori,K.PAPostCovariance,HObservationMatrix,K.RObservationNoiseCovariance,K.PHIStateTransitMatrix,K.QStateNoiseCovariance,1);

%=== Calculate Dilution Of Precision ======================================
%--- Initialize output ------------------------------------------------
dop     = zeros(1, 5);
A=HObservationMatrix(1:noSat,1:4);
%--- Calculate DOP ----------------------------------------------------
Q       = inv(A'*A);

dop(1)  = sqrt(trace(Q));                       % GDOP
dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
dop(4)  = sqrt(Q(3,3));                         % VDOP
dop(5)  = sqrt(Q(4,4));                         % TDOP


if nargout  == 4
    % RAIM parameters
    RAIMresult = RAIM(stateAPosteriori(1:4),KFMeasurement(1:noSat),HObservationMatrix(1:noSat,1:4),PAPrioriCovariance);
end
