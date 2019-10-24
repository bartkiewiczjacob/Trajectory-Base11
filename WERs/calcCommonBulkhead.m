function [wTot] = calcCommonBulkhead(id, p, mCH4, mLOX)

%% Description

% This function calculates the weight of common bulkhead tanks with given
% inner and outer diameters of LOX and CH4 tanks.


%% Initialization

ullage = .1 ; % recommended ullage percentage

rhoLOX = .0411; % density of LOX, lbm/in^3
rhoCH4 = .0150; % density of methane, lbm/in^3
rhoAl = .09754; % density of aluminum, lbm/in^3

wallCH4 = wallThickness(p, id); % minimum thickness of CH4 tank wall, in
wallLOX = wallThickness(p, id); % minimum thickness of LOX tank wall, in

odCH4 = id + wallCH4 * 2;
odLOX = id + wallLOX * 2;


%% Calculation

% recommended volumes of CH4 and LOX, in^3
vCH4 = mCH4 * (1 + ullage) / rhoCH4;
vLOX = mLOX * (1 + ullage) / rhoLOX;

% heights of tanks excluding domes, in
hCH4 = vCH4 ./ (pi * (id / 2) .^ 2);
hLOX = (vLOX - (4 / 3 * pi * (id / 2) .^ 3)) ./ (pi * (id / 2) .^ 2);

% weights of parts of tanks
domeCH4 = 2 / 3 * pi * rhoAl * ((odCH4 / 2) .^ 3 - (id / 2) .^ 3);
domeLOX = 4 / 3 * pi * rhoAl * ((odLOX / 2) .^ 3 - (id / 2) .^ 3);
wallCH4 = pi * .25 * rhoAl * (odCH4 .^ 2 - id .^ 2) .* hCH4;
wallLOX = pi * .25 * rhoAl * (odLOX .^ 2 - id .^ 2) .* hLOX;

wCH4 = domeCH4 + wallCH4; % weight of CH4 tank, lb
wLOX = domeLOX + wallLOX; % weight of LOX tank, lb

wTot = wCH4 + wLOX; % total weight of common bulkhead
