function [channels] = GPS_satpos(channels,activeChnList,transmitTime_increment)
%(transmitTime, prnList, eph) 
%SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
%given ephemeris EPH. Coordinates are calculated for each satellite in the
%list PRNLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph, settings);
%
%   Inputs:
%       transmitTime  - transmission time
%       prnList       - list of PRN-s to be processed
%       eph           - ephemeridies of satellites
%       settings      - receiver settings
%
%   Outputs:
%       satPositions  - positions of satellites (in ECEF system [X; Y; Z;])
%       satClkCorr    - correction of satellites clocks

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================
numOfSatellites = length(activeChnList);

% GPS constatns

gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate 
                                   % system

%--- Constants for satellite position calculation -------------------------
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986005e14;      % Earth's universal
                                   % gravitational parameter,
                                   % [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for ij = activeChnList
   eph = channels(ij).eph;
   transmitTime = transmitTime_increment(ij);
    
%% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime - eph.t_oc);

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(ij) = (eph.a_f2 * dt + eph.a_f1) * dt + ...
                         eph.a_f0 - ...
                         eph.T_GD;

    time = transmitTime - satClkCorr(ij);

%% Find satellite's position ----------------------------------------------

    %Restore semi-major axis
    A   = eph.sqrtA * eph.sqrtA;

    %Time correction
    tk  = check_t(time - eph.t_oe);

    %Initial mean motion
    n0  = sqrt(GM / A^3);
    %Mean motion
    n   = n0 + eph.deltan;

    %Mean anomaly
    M   = eph.M_0 + n * tk;
    %Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*gpsPi, 2*gpsPi);

    %Initial guess of eccentric anomaly
    E   = M;

    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph.e * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*gpsPi, 2*gpsPi);

    %Compute relativistic correction term
    dtr = F * eph.e * eph.sqrtA * sin(E);

    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph.e^2) * sin(E), cos(E)-eph.e);

    %Compute angle phi
    phi = nu + eph.omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*gpsPi);

    %Correct argument of latitude
    u = phi + ...
        eph.C_uc * cos(2*phi) + ...
        eph.C_us * sin(2*phi);
    %Correct radius
    r = A * (1 - eph.e*cos(E)) + ...
        eph.C_rc * cos(2*phi) + ...
        eph.C_rs * sin(2*phi);
    %Correct inclination
    i = eph.i_0 + eph.iDot * tk + ...
        eph.C_ic * cos(2*phi) + ...
        eph.C_is * sin(2*phi);

    %Compute the angle between the ascending node and the Greenwich meridian
    Omega = eph.omega_0 + (eph.omegaDot - Omegae_dot)*tk - ...
            Omegae_dot * eph.t_oe;
    %Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);

    %--- Compute satellite coordinates ------------------------------------
    channels(ij).satPositions(1) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
    channels(ij).satPositions(2) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
    channels(ij).satPositions(3) = sin(u)*r * sin(i);


      %--- Compute satellite velocities  %------------------------------------
    
     ekdot = n/(1.0 - eph.e*cos(E));

     tak = atan2( sqrt(1.0 - eph.e* eph.e)*sin(E), cos(E)- eph.e);
     takdot = sin(E)*ekdot*(1.0+eph.e*cos(tak))/(sin(tak)*(1.0-eph.e*cos(E)));


     
     ukdot = takdot +2.0*(eph.C_us*cos(2.0*u)-eph.C_uc*sin(2.0*u))*takdot;
     rkdot = A*eph.e*sin(E)*n/(1.0-eph.e*cos(E)) + 2.0*(eph.C_rs*cos(2.0*u)-eph.C_rc*sin(2.0*u))*takdot;
     ikdot = eph.iDot + (eph.C_is*cos(2.0*u)-eph.C_ic*sin(2.0*u))*2.0*takdot;

     xpk = r*cos(u);
     ypk = r*sin(u);

    xpkdot = rkdot*cos(u) - ypk*ukdot;
    ypkdot = rkdot*sin(u) + xpk*ukdot;

    omegak = Omega;
    omegakdot = (eph.omegaDot - Omegae_dot);

    xk = xpk*cos(omegak) - ypk*sin(omegak)*cos(i);
    yk = xpk*sin(omegak) + ypk*cos(omegak)*cos(i);
    zk =                   ypk*sin(i);

    xkdot = ( xpkdot-ypk*cos(i)*omegakdot )*cos(omegak)...
            - ( xpk*omegakdot+ypkdot*cos(i)-ypk*sin(i)*ikdot )*sin(omegak);
    ykdot = ( xpkdot-ypk*cos(i)*omegakdot )*sin(omegak)...
            + ( xpk*omegakdot+ypkdot*cos(i)-ypk*sin(i)*ikdot )*cos(omegak);
    zkdot = ypkdot*sin(i) + ypk*cos(i)*ikdot;

     
    channels(ij).satVelocity(1)= xkdot;
    channels(ij).satVelocity(2)= ykdot;
    channels(ij).satVelocity(3)= zkdot;
    channels(ij).satVelTOT     = sqrt(xkdot.^2+ykdot.^2+zkdot.^2)./1e3; %km/s
    
    
    
%% Include relativistic correction in clock correction --------------------
    channels(ij).satClkCorr = (eph.a_f2 * dt + eph.a_f1) * dt + ...
                         eph.a_f0 - ...
                         eph.T_GD + dtr;
%    if (~isempty(eph.A_0G))
%       channels(ij).GGTO = eph.A_0G + eph.A_1G*(eph.TOW-eph.t_0G+604800*mod((eph.WN-eph.WN_0G),64));
%    end
                     
end % for satNr = 1 : numOfSatellites
