%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
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
%
%Script initializes settings and environment of the software receiver.
%Then the processing is started.

%--------------------------------------------------------------------------
% CVS record:
% $Id: init.m,v 1.14.2.21 2006/08/22 13:46:00 dpl Exp $

%% Clean up the environment first =========================================
clearvars; close all; clc; fclose all;

format ('compact');
format ('long', 'g');

%--- Include folders with functions ---------------------------------------
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions
addpath ('RAIM',genpath('MPDD_functions'),'acq_track_functions','geoFunctions',genpath('include'),'Kalman Filter')%,'include for MPDD')                % RAIM functions

%% Print startup ==========================================================
fprintf(['\n',...
    'Welcome to:  softGNSS\n\n', ...
    'An open source GNSS SDR software project initiated by:\n\n', ...
    '              Danish GPS Center/Aalborg University\n\n', ...
    'The code was improved by GNSS Laboratory/University of Colorado.\n\n',...
    'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n',...
    'for details please read license details in the file license.txt. This\n',...
    'is free software, and  you  are  welcome  to  redistribute  it under\n',...
    'the terms described in the license.\n\n']);
fprintf('                   -------------------------------\n\n');

%% Initialize constants, settings =========================================
settings = initSettings();
settingsMPDD=InitsettingsMPDD(settings);

for datafileIndex=1:settings.N_of_datafile
    
    settings.datafileIndex=datafileIndex;
    if datafileIndex>1
        fclose all;
        clearvars -except datafileIndex settings settingsMPDD
    end
    %fileName=char(settings.fileName(datafileIndex));
    %fileName=char(settings.fileName{1,1}{1,datafileIndex});
    
    if settings.N_of_datafile==1 % Only 1 dataset or more?
        fileName=char(settings.fileName{1,1});
    else
        fileName=char(settings.fileName{1,1}{1,datafileIndex});      
    end
    
    
    % Bl=[0.5 1  2:2:20 30]
    % Bl=[0.4:0.1:0.6 1  2:2:20 25 30] analisi completa
    %for Bl=[2 6 10 20]%Bl=[0.5 1 2 6 8 10 20 30]
        
        %settings.dllNoiseBandwidth=Bl;
        
     %   for pseudorate=[500 1000] % 2000
            
            %settings.navSolPeriod=pseudorate;
            
            %% Generate plot of raw data and ask if ready to start processing =========
            try
                fprintf('Probing data (%s)...\n', fileName)
                probeData(fileName,settings);
            catch
                % There was an error, print it and exit
                errStruct = lasterror;
                disp(errStruct.message);
                disp('  (run setSettings or change settings in "initSettings.m" to reconfigure)')
                return;
            end
            
            disp('  Raw IF/IQ data plotted ')
            disp('  (run setSettings or change settings in "initSettings.m" to reconfigure)')
            disp(' ');
            gnssStart = 1; %input('Enter "1" to initiate GNSS processing or "0" to exit : ');
            navSolutions=[];
            
            if (gnssStart == 1)
                disp(' ');
                %start things rolling...
                postProcessing
                %%%%%%%%%%
                if ~isempty(navSolutions)
                    figure
                    plot(navSolutions.channel.rawP(:,2:end)','.-')
                    grid on
                    title('RAW \rho')
                    legend(num2str(navSolutions.channel.PRN(:,1)))
                    xlabel('Time [s]')
                    
                    %                 [firstSubFrame, activeChnList] = findPreambles(trackResults,settings);
                    %                 startPreamblems=firstSubFrame(1);
                    %                 discreteTime=[0:length(navSolutions.channel.rawP(end,:))-1];
                    %                 xtime=(discreteTime*settings.navSolPeriod+startPreamblems)/1e3;
                    %                 plotNavigation2(navSolutions,settings,xtime);
                    
                    
                end
            end
            
            
%             close all
%             clearvars -except settings inddataset pseudorate Bl
%             clc
            
        %end
        
 %   end
    
end