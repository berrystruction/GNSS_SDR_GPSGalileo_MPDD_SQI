function plotNavigation(navSolutions,settings,myfolder)
%Functions plots variations of coordinates over time and a 3D position
%plot. It plots receiver coordinates in UTM system or coordinate offsets if
%the true UTM receiver coordinates are provided.
%
%plotNavigation(navSolutions, settings)
%
%   Inputs:
%       navSolutions    - Results from navigation solution function. It
%                       contains measured pseudoranges and receiver
%                       coordinates.
%       settings        - Receiver settings. The true receiver coordinates
%                       are contained in this structure.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

% CVS record:
% $Id: plotNavigation.m,v 1.1.2.25 2006/08/09 17:20:11 dpl Exp $

%% Plot results in the necessary data exists ==============================
if (~isempty(navSolutions))
    
    %% If reference position is not provided, then set reference position
    %% to the average postion
    if isnan(settings.truePosition.E) || isnan(settings.truePosition.N) ...
            || isnan(settings.truePosition.U)
        
        %=== Compute mean values ==========================================
        % Remove NaN-s or the output of the function MEAN will be NaN.
        refCoord.E = mean(navSolutions.E(~isnan(navSolutions.E)));
        refCoord.N = mean(navSolutions.N(~isnan(navSolutions.N)));
        refCoord.U = mean(navSolutions.U(~isnan(navSolutions.U)));
        
        %Also convert geodetic coordinates to deg:min:sec vector format
        %         meanLongitude = dms2mat(degrees2dms(...
        %             mean(navSolutions.longitude(~isnan(navSolutions.longitude)))), -5);
        %         meanLatitude  = dms2mat(degrees2dms(...
        %             mean(navSolutions.latitude(~isnan(navSolutions.latitude)))), -5);
        meanLongitude = degrees2dms(...
            mean(navSolutions.longitude(~isnan(navSolutions.longitude))));
        meanLatitude  = degrees2dms(...
            mean(navSolutions.latitude(~isnan(navSolutions.latitude))));
        
        refPointLgText = ['Mean Position\newline  Lat: ', ...
            num2str(meanLatitude(1)), '{\circ}', ...
            num2str(meanLatitude(2)), '{\prime}', ...
            num2str(meanLatitude(3)), '{\prime}{\prime}', ...
            '\newline Lng: ', ...
            num2str(meanLongitude(1)), '{\circ}', ...
            num2str(meanLongitude(2)), '{\prime}', ...
            num2str(meanLongitude(3)), '{\prime}{\prime}', ...
            '\newline Hgt: ', ...
            num2str(mean(navSolutions.height(~isnan(navSolutions.height))), '%+6.1f')];
    else
        refPointLgText = 'Reference Position';
        refCoord.E = settings.truePosition.E;
        refCoord.N = settings.truePosition.N;
        refCoord.U = settings.truePosition.U;
    end
    
    figureNumber = round(50+rand()*150);
    % The 300 is chosen for more convenient handling of the open
    % figure windows, when many figures are closed and reopened. Figures
    % drawn or opened by the user, will not be "overwritten" by this
    % function if the auto numbering is not used.
    
    %=== Select (or create) and clear the figure ==========================
    figure(figureNumber);
    clf   (figureNumber);
    set   (figureNumber, 'Name', 'Navigation solutions');
    
    %--- Draw axes --------------------------------------------------------
    handles(1, 1) = subplot(4, 2, 1 : 4);
    handles(3, 1) = subplot(4, 2, [5, 7]);
    handles(3, 2) = subplot(4, 2, [6, 8]);
    
    %% Plot all figures =======================================================
    
    %--- Coordinate differences in UTM system -----------------------------
    xtime=([0:length(navSolutions.E(end,:))-1]+floor(settings.msToProcess/settings.navSolPeriod-length(navSolutions.E(end,:))))*settings.navSolPeriod/1e3;
    plot(handles(1, 1), xtime, (navSolutions.E - refCoord.E),'.-', ...
        xtime, (navSolutions.N - refCoord.N),'.-',...
        xtime, (navSolutions.U - refCoord.U),'.-');
    
    meansquareErr=sqrt([mean((navSolutions.E(2:end-1)-refCoord.E).^2), mean((navSolutions.N(2:end-1)-refCoord.N).^2), mean((navSolutions.U(2:end-1)-refCoord.U).^2)]);
    avgDist=mean(sqrt((navSolutions.E(2:end-1)-refCoord.E).^2+(navSolutions.N(2:end-1)-refCoord.N).^2+(navSolutions.U(2:end-1)-refCoord.U).^2));

    disp('ENU errors and MSE:')
    disp([num2str([meansquareErr avgDist])])
    
    title (handles(1, 1), ['Coordinates variations in UTM system. Avg. ENU Distance ' num2str(avgDist) ' m. MSE for E=' num2str(meansquareErr(1)) ' m, MSE for N=' num2str(meansquareErr(2)) ' m, MSE for U=' num2str(meansquareErr(3)) ' m']);
    legend(handles(1, 1), 'E', 'N', 'U');
    xlabel(handles(1, 1), ['Time [s] (Measurement period: ', ...
        num2str(settings.navSolPeriod), ' ms)']);
    ylabel(handles(1, 1), 'Variations (m)');
    grid  (handles(1, 1));
    %     axis  (handles(1, 1), 'tight');
    allnavSol=[navSolutions.E-refCoord.E, navSolutions.N-refCoord.N, navSolutions.U-refCoord.U];
    axis  (handles(1, 1), [0 settings.msToProcess/1e3 min(allnavSol) max(allnavSol)]);
    
    %--- Position plot in UTM system --------------------------------------
    plot3 (handles(3, 1), navSolutions.E - refCoord.E, ...
        navSolutions.N - refCoord.N, ...
        navSolutions.U - refCoord.U, '+');
    hold  (handles(3, 1), 'on');
    %Plot the reference point
    plot3 (handles(3, 1), 0, 0, 0, 'r+', 'LineWidth', 1.5, 'MarkerSize', 10);
    hold  (handles(3, 1), 'off');
    
    view  (handles(3, 1), 0, 90);
    axis  (handles(3, 1), 'equal');
    grid  (handles(3, 1), 'minor');
    
    legend(handles(3, 1), 'Measurements', refPointLgText);
    
    title (handles(3, 1), 'Positions in UTM system (3D plot)');
    xlabel(handles(3, 1), 'East (m)');
    ylabel(handles(3, 1), 'North (m)');
    zlabel(handles(3, 1), 'Upping (m)');
    
    %--- Satellite sky plot -----------------------------------------------
    skyPlot(handles(3, 2), ...
        navSolutions.channel.az, ...
        navSolutions.channel.el, ...
        navSolutions.channel.PRN(:, 1));
    
    title (handles(3, 2), ['Sky plot (mean PDOP: ', ...
        num2str(mean(navSolutions.DOP(2,:))), ')']);
    %% plot of GDOP in time    
    if isfield(navSolutions,'DOP')==1
        title (handles(3, 2), ['Sky plot (mean PDOP: ', ...
        num2str(mean(navSolutions.DOP(2,:))), ')']);
    
        figure
        if navSolutions.TypeofSol==0 % if =0, WLS solution
            subplot(221)
            plot(xtime,navSolutions.DOP(1,:),xtime,navSolutions.WDOP(1,:),'.-'), grid on, xlabel('Time [s]'), title('GDOP trend'), legend('GDOP','WGDOP')
            subplot(222)
            plot(xtime,navSolutions.DOP(2,:),xtime,navSolutions.WDOP(2,:),'.-'), grid on, xlabel('Time [s]'), title('PDOP vs \sigma_V_W'), legend('PDOP','WPDOP')
            subplot(223)
            plot(xtime,navSolutions.DOP(3,:),xtime,navSolutions.WDOP(3,:),'.-'), grid on, xlabel('Time [s]'), title('HDOP vs \sigma_H_W'), legend('HDOP','WHDOP')
            subplot(224)
            plot(xtime,navSolutions.DOP(4,:),xtime,navSolutions.WDOP(4,:),'.-'), grid on, xlabel('Time [s]'), title('VDOP vs \sigma_V_W'), legend('VDOP','WVDOP')
        else
            plot(xtime,navSolutions.DOP(1,:),'.-'), grid on, xlabel('Time [s]'), title('GDOP trend')
        end
    else
        title('Sky Plot');
    end
    
    %% Protection Level vs Time plot
    if settings.RAIM.enableRAIM==true
        figure
        subplot(221)
        plot(xtime,[navSolutions.PL(:).HPL],'r',xtime,[navSolutions.PL(:).HPL2],'bo'), title('Horizontal Protection Level'), grid on, xlabel('Time [s]'), ylabel('HPL [m]')
        subplot(222)
        plot(xtime,[navSolutions.PL(:).VPL],'r',xtime,[navSolutions.PL(:).VPL2],'bo'), title('Vertical Protection Level'), grid on, xlabel('Time [s]'), ylabel('VPL [m]')
        subplot(223)
        plot(xtime,[navSolutions.PL(:).HRMS],'r',xtime,[navSolutions.PL(:).HRMS2],'b'), title('Horizontal Positioning accuracy'), grid on, xlabel('Time [s]'), ylabel('HRMS [m]')
        subplot(224)
        plot(xtime,[navSolutions.PL(:).VRMS],'r',xtime,[navSolutions.PL(:).VRMS2],'b'), title('Vertical Positioning accuracy'), grid on, xlabel('Time [s]'), ylabel('VRMS [m]')
    end
    %% Writing on .kml file to visualize the track on Earth 
    KML_googleEarth([myfolder '.kml'],navSolutions.latitude,navSolutions.longitude,navSolutions.height)
    
else
    disp('plotNavigation: No navigation data to plot.');
end % if (~isempty(navSolutions))
fprintf('\n');

