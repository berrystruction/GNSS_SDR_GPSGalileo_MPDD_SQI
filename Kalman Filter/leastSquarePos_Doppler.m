function [pos_fin, el, az, dop, RAIMresult] = leastSquarePos_Doppler(satpos, obs,channels,doppler, settings)
%Function calculates the Least Square Solution.
%
%[pos, el, az, dop] = leastSquarePos(satpos, obs, channels, doppler);
%
%   Inputs:
%       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)
%       obs         - Observations - the pseudorange measurements to each
%                   satellite:
%                   (e.g. [20000000 21000000 .... .... .... .... ....])
%       channels    - receiver settings
%       doppler     - doppler freq
%
%   Outputs:
%       pos_fin     - receiver position, velocy and receiver clock error and clock drift
%                   (in ECEF system: [X, Y, Z, Clk_T, dX, dY, dZ, dClk_T ])
%       el          - Satellites elevation angles (degrees)
%       az          - Satellites azimuth angles (degrees)
%       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: leastSquarePos.m,v 1.1.2.12 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

%=== Initialization =======================================================
nmbOfIterations = 7;

f0_L1 = settings.f0_L1; % portante L1 in Hz

dtr     = pi/180;
pos     = zeros(4, 1);
pos_fin = zeros(8, 1);
X       = satpos;
nmbOfSatellites = size(satpos, 2);
Rot_X = zeros(3,nmbOfSatellites);

A       = zeros(nmbOfSatellites, 4);
omc     = zeros(nmbOfSatellites, 1);
az      = zeros(1, nmbOfSatellites);
el      = az;
sol = settings.c;
%=== Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations
    
    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X(:,i) = X(:, i);
            trop = 2;
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / sol ;
            
            %--- Correct satellite position (do to earth rotation) --------
            Rot_X(:,i) = e_r_corr(traveltime, X(:, i));
            %--- Find the elevation angle of the satellite ----------------
            [az(i), el(i), dist] = topocent(pos(1:3, :), Rot_X(:,i) - pos(1:3, :));
            useTropCorr = 1;
            if (useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                    0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
                %file_1 = fopen('output1.txt','a');
                %fprintf(file_1, '%f\n', trop);
                %fclose(file_1)
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end
        end % if iter == 1 ... ... else
        
        %--- Apply the corrections ----------------------------------------
        omc(i) = (obs(i) - norm(Rot_X(:,i) - pos(1:3), 'fro') - pos(4) - trop);
        
        %--- Construct the A matrix ---------------------------------------
        A(i, :) =  [ (-(Rot_X(1,i) - pos(1))) / obs(i) ...
            (-(Rot_X(2,i) - pos(2))) / obs(i) ...
            (-(Rot_X(3,i) - pos(3))) / obs(i) ...
            1 ];
    end % for i = 1:nmbOfSatellites
    
    % These lines allow the code to exit gracefully in case of any errors
    if rank(A) ~= 4
        pos     = zeros(1, 4);
        dop     = zeros(1, 5);
        return
    end
    
    %--- Find position update ---------------------------------------------
    x   = A \ omc;
    %keyboard
    %--- Apply position update
    %--------------------------------------------
    pos = pos + x;
    
end % for iter = 1:nmbOfIterations

pos = pos';

%=== Calculate velocity estimation ========================================
%=== "Understanding GPS - Kaplan" ppg. 48-51
%=== Defining column vectors used in computation ==========================
stimatore =2;
if stimatore == 1 % derivata
    user_vel = (pos - pos_old)/settings.navSolPeriod;
elseif stimatore == 2 % metodo doppler
    % Ricostruisco A
    A_corretta = zeros(size(A));
    for i = 1:nmbOfSatellites
        distVec = Rot_X(:,i)' - pos(1:3);
        distanza = sqrt(distVec * distVec');
        A_corretta(i, :) =  [distVec/distanza, 1];
    end
    
    % 	if ~isnan(sum(X_old(:))) || velEst == 2 % X_old ??? NaN solo alla prima iterazione
    % 		% TODO settings.navSolPeriod ??? corretto?
    %
    sat_radial_vel_and_doppler = zeros(nmbOfSatellites,1);
    for i = 1:nmbOfSatellites
        % 			if velEst == 1
        % 			% sat_vel ??? una colonna, A(i,1:3) ??? una riga -> prodotto scalare
        % 			sat_vel = (Rot_X(:,i) - X_old(:,i))/(settings.navSolPeriod/1000); % m/s
        % 			elseif velEst == 2
        sat_vel = channels(i).satVelocity(1:3)';
        % 			else
        % 				fprintf('\nStimatore della velocit??? del satellite errato!!!\n');
        % 				return
        % 			end
        sat_radial_vel_and_doppler(i) = sol / f0_L1 * doppler(i) +...
            A_corretta(i,1:3) * sat_vel; % sul libro ??? il vettore d
        
    end
    % risoluzione del sistema lineare
    user_vel = A_corretta \ sat_radial_vel_and_doppler; % contiene anche il drift rate del clock del rx
    
    
    pos_fin(1:4) = pos;
    pos_fin(5:8) = user_vel;
    % 	else
    % 		user_vel = NaN(4,1); % alla prima iterazione non stimo la velocit???
    % 	end
else
    fprintf('Metodo inesistente\n');
    return;
end

%=== Calculate Dilution Of Precision ======================================
% if nargout  == 4
%--- Initialize output ------------------------------------------------
dop     = zeros(1, 5);

%--- Calculate DOP ----------------------------------------------------
Q       = inv(A'*A);

dop(1)  = sqrt(trace(Q));                       % GDOP
dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
dop(4)  = sqrt(Q(3,3));                         % VDOP
dop(5)  = sqrt(Q(4,4));                         % TDOP

%=== RAIM parameters======================================
if nargout  == 5
    RAIMresult = RAIM(x,omc,A,nmbOfSatellites);   
end
