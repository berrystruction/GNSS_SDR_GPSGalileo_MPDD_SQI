function plotNavigationNoref(navSolutions1,navSolutions2)
%Functions plots variations of coordinates over time and a 3D position
%plot.
%plotNavigation(navSolutions, settings)
%
%   Inputs:
%       navSolutions    - Results from navigation solution function. It
%                       contains measured pseudoranges and receiver
%                       coordinates.
%       settings        - Receiver settings. The true receiver coordinates
%                       are contained in this structure.

%--------------------------------------------------------------------------


%% Plot results in the necessary data exists ==============================
if (~isempty(navSolutions1))
    
    figureNumber = round(50+rand()*150);
    % The 300 is chosen for more convenient handling of the open
    % figure windows, when many figures are closed and reopened. Figures
    % drawn or opened by the user, will not be "overwritten" by this
    % function if the auto numbering is not used.
    
    %=== Select (or create) and clear the figure ==========================
    figure(figureNumber);
    clf   (figureNumber);
    set   (figureNumber, 'Name', 'ENU Navigation solutions');
    if nargin==2
        
        %--- Position plot in UTM system --------------------------------------
        plot3 ( navSolutions1.E-navSolutions2.E, ...
            navSolutions1.N-navSolutions2.N, ...
            navSolutions1.U-navSolutions2.U, '.');
        %hold  ( 'on');
    else
        %--- Position plot in UTM system --------------------------------------
        plot3 (navSolutions1.E, ...
            navSolutions1.N, ...
            navSolutions1.U, 'o');
    end
    view  ( 0, 90);
    axis  ( 'equal');
    grid  ( 'minor');
    
    title ( 'Positions in UTM system (3D plot)');
    xlabel( 'East (m)');
    ylabel( 'North (m)');
    zlabel( 'Upping (m)');
    
    
    
else
    disp('plotNavigation: No navigation data to plot.');
end % if (~isempty(navSolutions))
fprintf('\n');

