function [pos, el, az, dop, wdop, RAIMresult] = WLS(satpos, obs, settings, W, eph_for_IONO, TOW)
%Function calculates the Least Square Solution.
%
%[pos, el, az, dop] = leastSquarePos(satpos, obs, settings);
%
%   Inputs:
%       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)
%       obs         - Observations - the pseudorange measurements to each
%                   satellite:
%                   (e.g. [20000000 21000000 .... .... .... .... ....])
%       settings    - receiver settings
%
%   Outputs:
%       pos         - receiver position and receiver clock error
%                   (in ECEF system: [X, Y, Z, dt])
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

%invW=W\eye(size(W)); % inverse of covariance matrix
invW=inv(W);

dtr     = pi/180;
pos     = zeros(4, 1);
X       = satpos;
nmbOfSatellites = size(satpos, 2);

H       = zeros(nmbOfSatellites, 4);
omc     = zeros(nmbOfSatellites, 1);
az      = zeros(1, nmbOfSatellites);
el      = az;

if eph_for_IONO.iflag==1 % if Ionospheric parameters exist
    Alpha=[eph_for_IONO.alpha_0 eph_for_IONO.alpha_1 eph_for_IONO.alpha_2 eph_for_IONO.alpha_3]';
    Beta=[eph_for_IONO.beta_0 eph_for_IONO.beta_1 eph_for_IONO.beta_2 eph_for_IONO.beta_3]';
end

%=== Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations
    
    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);
            trop = 2;
            IONOcorr=0;            
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c ;
            
            %--- Correct satellite position (do to earth rotation) --------
            Rot_X = e_r_corr(traveltime, X(:, i));
            
            %--- Find the elevation angle of the satellite ----------------
            [az(i), el(i), dist] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));
            
            if (settings.useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                    0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end
            
            % IONO corrections
            if eph_for_IONO.iflag==1 % if Ionospheric params exist, compute the correction
                IONOcorr=Error_Ionospheric_Klobuchar(pos(1:3)',Rot_X',Alpha,Beta,TOW);
                IONOcorr=IONOcorr*settings.c;
            else
                IONOcorr=0;
            end
            
        end % if iter == 1 ... ... else
        
        %--- Apply the corrections ----------------------------------------
        omc(i) = (obs(i) - norm(Rot_X - pos(1:3), 'fro') - pos(4) - trop - IONOcorr);
        
        %--- Construct the A matrix ---------------------------------------
        H(i, :) =  [ (-(Rot_X(1) - pos(1))) / obs(i) ...
            (-(Rot_X(2) - pos(2))) / obs(i) ...
            (-(Rot_X(3) - pos(3))) / obs(i) ...
            1 ];
    end % for i = 1:nmbOfSatellites
    
    % These lines allow the code to exit gracefully in case of any errors
    if rank(H) ~= 4
        pos     = zeros(1, 4);
        return
    end
    
    %--- Find position update/WLS solution  ---------------------------------------------
    x   = (H'*invW*H)^(-1)*H'*invW*omc;
    %--- Apply position update --------------------------------------------
    pos = pos + x;
    
    %     if nmbOfIterations==iter
    %         omcc=omc;
    %     end
    
end % for iter = 1:nmbOfIterations

% RAIM extension
% RAIMresult.pos=pos;
% RAIMresult.x=x;
% RAIMresult.omc=omc;
% RAIMresult.nmbOfSatellites=nmbOfSatellites;
% RAIMresult.A=A;
% RAIMresult.chlist=0;

if settings.RAIM.enableRAIM==true
    RAIMresult = RAIM(x,omc,H,nmbOfSatellites);
else
    RAIMresult=[];
end
%%

pos = pos';

%=== Calculate Dilution Of Precision ======================================
if nargout  == 6
    %--- Initialize output ------------------------------------------------
    dop     = zeros(1, 5);
    
    %--- Calculate DOP ----------------------------------------------------
    Q       = inv(H'*H);
    
    dop(1)  = sqrt(trace(Q));                       % GDOP
    dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
    dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
    dop(4)  = sqrt(Q(3,3));                         % VDOP
    dop(5)  = sqrt(Q(4,4));                         % TDOP
    
    
    %--- Initialize output ------------------------------------------------
    wdop     = zeros(1, 5);
    %--- Calculate DOP ----------------------------------------------------
    Q       = inv(H'*invW*H);
    
    wdop(1)  = sqrt(trace(Q));                       % WGDOP
    wdop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % WPDOP
    wdop(3)  = sqrt(Q(1,1) + Q(2,2));                % WHDOP
    wdop(4)  = sqrt(Q(3,3));                         % WVDOP
    wdop(5)  = sqrt(Q(4,4));                         % WTDOP
end

