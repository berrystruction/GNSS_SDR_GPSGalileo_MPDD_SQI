function [trueRange,el,az] = pseudo2true(pdr, satClkCorr, pos, dt, satpos,c)

%c = 299792458;
noSat = size(satpos,2);

[ert(1) ert(2) ert(3)] = xyz2llh(pos(1),pos(2),pos(3));
h = ert(3);

trop = zeros(noSat,1);
az = trop;
el = trop;

% Calcolo dell'elevazione dei satelliti
for ind = 1:noSat
    [az(ind),el(ind)] = ComputeAzEl(satpos(:,ind),pos);
    [trop(ind), Sigma2Tropo] = TropoHopfieldModel(el(ind),h);
end

trueRange = pdr + satClkCorr.' * c - dt + trop;