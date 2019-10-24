function [wTot] = calcCoAx(idLOX, p, mCH4, mLOX)

%% Description

% This function calculates the weight of Co-Ax tanks with given inner and
% outer diameters of LOX and CH4 tanks.


%% Initialization

ullage = .1 ; % recommended ullage percentage
dP = 50; % pressure differential for shared wall

rhoLOX = .0411; % density of LOX, lbm/in^3
rhoCH4 = .0150; % density of methane, lbm/in^3
rhoAl = .09754; % density of aluminum, lbm/in^3


%% Calculation

% recommended volumes of CH4 and LOX, in^3
vCH4 = mCH4 * (1 + ullage) / rhoCH4;
vLOX = mLOX * (1 + ullage) / rhoLOX;
vTot = vCH4 + vLOX; % total volume of propellants, in^3

% wall thickness, in
tLOX = wallThickness(p, idLOX); % minimum thickness of LOX tank wall, in
odLOX = idLOX + 2 * tLOX;

% height of tanks excluding domes, in
odCH4 = .5 * idLOX;
approx = idLOX;
error = 0.1;
while max(abs(odCH4 - approx)) > error
    odCH4 = (odCH4 + approx) / 2;
    h = abs((vTot - (4 / 3 * pi * (idLOX / 2) .^ 2 .* (odCH4 / 2))) ./ (pi * (idLOX / 2) .^ 2));
    tCH4 = real(wallThickness(dP, odCH4)); % thickness approximately the same when OD is used for ID
    idCH4 = real(odCH4 - 2 * tCH4);
    approx = ((vCH4 - pi * (idCH4 / 2) .^ 2 .* h) / (4 / 3 * pi)) .^ (1 / 3) + 2 * tCH4;
end
odCH4 = approx;
idCH4 = odCH4 - 2 * tCH4;

% weights of parts of tanks
domeCH4 = 4 / 3 * pi * rhoAl * ((odCH4 / 2) .^ 3 - (idCH4 / 2) .^ 3);
domeLOX = 4 / 3 * pi * rhoAl * ((odLOX / 2) .^ 2 .* ((odLOX - idLOX + odCH4) / 2) - (idLOX / 2) .^ 2 .* (odCH4 / 2));
wallCH4 = pi * .25 * rhoAl * (odCH4 .^ 2 - idCH4 .^ 2) .* h;
wallLOX = pi * .25 * rhoAl * (odLOX .^ 2 - idLOX .^ 2) .* h;

wCH4 = domeCH4 + wallCH4; % weight of CH4 tank, lb
wLOX = domeLOX + wallLOX; % weight of LOX tank, lb

wTot = wCH4 + wLOX; % total weight of common bulkhead
