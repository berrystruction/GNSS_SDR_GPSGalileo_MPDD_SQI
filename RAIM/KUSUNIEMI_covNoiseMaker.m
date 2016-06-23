function [covarianceMTX, covariancePosVelocityMTX]  = covNoiseMaker(trackres,currMeasNr,lenChList)
%COVNOISEMAKER
% Heidi Kuusniemi model
% Thesis: User-Level Reliability and Quality Monitoring in Satellite-Based Personal
% Navigation

CN0sorted=zeros(1,lenChList);
for index=1:lenChList
    CN0sorted(index)=trackres(index).CN0.CNo_SNV(currMeasNr);
end
%cov noise matrix for NO MP dataset 5 min different CN0
%load covNoiseMtx.mat

%cov noise matrix for MP dataset 20 min
Nmeasurement=length(CN0sorted);
% CN0sorted=45*ones(1,Nmeasurement);

% Lightly degreded environment paramenters
a=10;  % [m^2]
b=150^2; % [m^2*Hz]
% Heavily degreded environment paramenters
% a=500;  % [m^2]
% b=10^6; % [m^2*Hz]
sigma2_i=a+b*10.^(-CN0sorted/10);
covarianceMTX=zeros(Nmeasurement);
aus=eye(Nmeasurement);
for index=1:Nmeasurement
    covarianceMTX(index,index)=sigma2_i*aus(:,index);
end

if nargout==2
    % Lightly degreded environment paramenters
    a=0.01;  % [m^2/s^2]
    b=25; % [m^2/s^2*Hz]
    % Heavily degreded environment paramenters
%     a=0.001;  % [m^2/s^2]
%     b=40^6; % [m^2/s^2*Hz]
    sigma2_i=a+b*10.^(-CN0sorted/10);
    covariancePosVelocityMTX=zeros(2*Nmeasurement);
    covariancePosVelocityMTX(1:Nmeasurement,1:Nmeasurement)=covarianceMTX;
    for index=Nmeasurement+1:2*Nmeasurement
        covariancePosVelocityMTX(index,index)=sigma2_i*aus(:,index-Nmeasurement);
    end     
end

% Metodo di Davide- frobenius norm
%sigmaeff=sqrt(trace((covw*covw'))/(Nmeasurement-4));
%sigma2=31.9410;
% sigma2=8.73169075152966;
% covarianceMTX=sigma2*eye(Nmeasurement);

% tabella=[1  4  6  11 16 18 19 21 22 27;...
%          42 45 45 46 44 43 46 40 43 50];

end

