function K = initKalmanGPS

% ---- Kalman data --------------------------------------------------------
% Nel 3D ho i seguenti stati:
% x y z t dx dy dz dt
K.Nstate = 8;
K.stateAPosteriori = zeros(K.Nstate,1);

% posizione e velocità (spazio e tempo)
%K.posSigma2 = 10e01; % m^2, varianza delle misure
%K.velSigma2 = 1e-2; % m^2/s^2, varianza delle misure

K.PAPostCovariance = 5e2*diag([2e0 2e0 2e0 2e6 2e0 2e0 2e0 2e1]); % è l'indice di affidabilità delle stime 5e0

K.PHIStateTransitMatrix = eye(8);
K.PHIStateTransitMatrix(1:4,5:8) = eye(4);

Sp = 5e-2; Sf = 1e-1; Sg = 1e-2;

%K.QStateNoiseCovariance = 1e2*diag([10*Sp/3 10*Sp/3 10*Sp/3 (Sf+Sg/2) 1*Sp 1*Sp 1*Sp Sg]); % !ORIGINAL LINE! 1e0
K.QStateNoiseCovariance = 1e2*diag([10*Sp/3 10*Sp/3 10*Sp/3 (Sf+Sg/2) 1*Sp 1*Sp 1*Sp Sg]);

for i = 0:2
    K.QStateNoiseCovariance(5+i,1+i) = Sp/2;
    K.QStateNoiseCovariance(1+i,5+i) = Sp/2;
end
clear i;
K.QStateNoiseCovariance(8,4) = Sg/2;
K.QStateNoiseCovariance(4,8) = Sg/2;

% modello nominale
K.velocity = [0; 0; 0];

% K.fattorePOS = 1.8;
% K.fattoreVEL = 1.5;
% K.soglia = 38;

K.h = 330;
K.hvar = 14;