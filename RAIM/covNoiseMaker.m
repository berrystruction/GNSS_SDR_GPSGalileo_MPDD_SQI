function [covarianceMTX, covariancePosVelocityMTX, cumulativepseudo]  = covNoiseMaker(trackres,currMeasNr,lenChList,pseudoNew,pseudocumul,SQI,SQIthreshold,SQIpenalty)
%COVNOISEMAKER
% Heidi Kuusniemi model and mine
% Thesis: User-Level Reliability and Quality Monitoring in Satellite-Based Personal
% Navigation


% Kuusniemi parameters
CN0sorted=zeros(1,lenChList);
Nmeasurement=lenChList;
for index=1:lenChList
    CN0sorted(index)=trackres(index).CN0.CNo_SNV(currMeasNr);
end
% Lightly degraded environment paramenters
a=10;  % [m^2]
b=150^2; % [m^2*Hz]
% Heavily degreded environment paramenters
% a=500;  % [m^2]
% b=10^6; % [m^2*Hz]
sigma2_i=a+b*10.^(-CN0sorted/10);
covarianceMTX=zeros(Nmeasurement);
aus=eye(Nmeasurement);

for index=1:lenChList
    
    if SQI(index)>SQIthreshold && currMeasNr>10 && 1==0
        
        % DA CONTINUARE...
        cumulativepseudo=(pseudocumul+pseudoNew)/M;
        
    else
        
        %for index=1:Nmeasurement
        covarianceMTX(index,index)=sigma2_i*aus(:,index);
        %end
    end
    
end




% Creo la matrice di pesi con SQI
if 1==1 %&& SQIpenalty>0%==true % SQI Penalty in covariance noise matrix
    W=eye(lenChList);
    % % % %                 Ei=satElev(activeChnList)'*pi/180;
    % % % %                 %wsummation=sum(SQIchannels([trackResults(activeChnList).PRN],currMeasNr)+(cos(pi/2-Ei).^2+sin(pi/2-Ei).^2));
    % % % %                 wsummation=cos(pi/2-Ei).^2+0.3*sin(pi/2-Ei).^2;
    for index=1:lenChList
        % % % %                     %Ei= satElev(activeChnList(index))*pi/180;% Ei=0;
        % % % %                     %W(index,index)=(SQIchannels([trackResults(activeChnList(index)).PRN],currMeasNr)+(cos(pi/2-Ei).^2+sin(pi/2-Ei).^2))/wsummation*W(index,index);
        % % % %                     %W(index,index)=wsummation(index)*W(index,index);
        
        W(index,index)=SQIpenalty*(1-SQI(index))*covarianceMTX(index,index)+covarianceMTX(index,index);
    end
    covarianceMTX=W;
end








if nargout>1%nargout==2
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


cumulativepseudo=0;
end

