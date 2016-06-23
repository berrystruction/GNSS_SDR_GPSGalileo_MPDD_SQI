function navSolComparison(calcpositions,RefPositions,navSolPeriod,calcpositions2)
% Comparison between calculated positions and reference ones by aligning on
%correct TOW
% INPUT:
% - RefPositions: reference positions in ECEF coordinates and
% information about TOW
% - calcpositions: calculated posititons in ECEF and TOW information

figure
[~,maxpos]=max(abs(calcpositions.TOW(1)-RefPositions(:,4))<1e-5); % aligning of navigation solutions
[~,endpos]=max(abs(calcpositions.TOW(end)-RefPositions(:,4))<1e-5); % aligning of navigation solutions

%lenmax=min(length(calcpositions.X),length(RefPositions))+1;

% RefPositionsllh=zeros(lenmax);
% [RefPositionsllh(:,1),RefPositionsllh(:,2),RefPositionsllh(:,3)]=cart2geo(RefPositions(maxpos:navSolPeriod:lenmax*1000,1),...
%     RefPositions(maxpos:navSolPeriod:lenmax*1000,2),...
%     RefPositions(maxpos:navSolPeriod:lenmax*1000,3),5);
% utmzone=findUtmZone(RefPositionsllh(1,1),RefPositionsllh(1,2));
% [E,N,U]=cart2utm(RefPositionsllh(:,1),RefPositionsllh(:,2),RefPositionsllh(:,3),utmzone); % Reference positions

refX=RefPositions(maxpos:navSolPeriod:endpos,1);
refY=RefPositions(maxpos:navSolPeriod:endpos,2);
refZ=RefPositions(maxpos:navSolPeriod:endpos,3);

subplot(211)
plot(refX-calcpositions.X'), hold on
plot(refY-calcpositions.Y'), hold on
plot(refZ-calcpositions.Z')

grid minor
title('Reference positions Vs Calculated positions')
xlabel('Time')
legend('X','Y','Z')
subplot(212)
plot(RefPositions(maxpos:navSolPeriod:endpos,4)-calcpositions.TOW')
grid on
title('Aligned TOW')


if nargin==4
    figure,
    plot3(calcpositions2.X'-refX,calcpositions2.Y'-refY,calcpositions2.Z'-refZ,'.'), hold on
    plot3(calcpositions.X'-refX,calcpositions.Y'-refY,calcpositions.Z'-refZ,'o')
    title('LS vs KF Navigation solution')
    xlabel('East [m]'), grid on, ylabel('North [m]'), zlabel('Up [m]'), legend('LS','KF')
end

end

