
function [mass,batt] = BatteryModel(Vdc,Ebatt,Idc,D,tank,intertank)

%%
%
% Call: [mass,batt] = BatteryModel(Vdc,Ebatt,Idc,D);
%
% Input:
%
%  Vdc          = DC bus voltage (rated voltage) [V]
%  Ebatt        = Battery energy [Wh]
%  Idc          = DC bus current [A]
%  D            = Design specs parameters structure
%   D.name      = Battery Name
%   D.Vdc       = Battery Voltage per cell [V]
%   D.Idc       = Battery DC Current cell [A]
%   D.C         = Battery Charge per cell [Ah]
%   D.Mass      = Battery Mass per cell [kg]
%
% Output:
%
%  mass         = Battery mass [lbm]
%  batt         = Battery parameters output structure
%   batt.Name   = Name of the battery
%   batt.Ns     = Number of cells in series
%   batt.Np     = Number of cells in parallel


%% Pull Variables from the Struture
Dname = extractfield(D,'name'); % battery cell name
DVdc = extractfield(D,'Vdc'); % battery cell voltage [Vdc]
DIdc = extractfield(D,'Idc'); % battery cell current {Idc]
DC = extractfield(D,'C'); % battery cell charge [Ah]
Dmass = extractfield(D,'mass')*2.20462; % battery cell mass [lbm]

%% Size up Ebatt to account for voltage dropoff

Ebatt = 1.15*Ebatt;

%% Number of batteries in each Battery Pack
Ns = ceil(Vdc./DVdc); % cells in series
Np = ceil(Idc./DIdc); % cells in parallel

%% Add batteries in Series to Account for voltage drop

K =  21.2; % Resisitivity of aluminum [ohms-cmil/ft]

% NEC ampacity ratings (see:
% https://www.omnicable.com/docs/default-source/technical-resources/nec-ampacity-data.pdf?sfvrsn=2)
currentrat = [292 328 364 395 458 514 570]; 
gauge = [250 300 350 400 500 600 700];

% selects rating of the cable
rating = find(currentrat > Idc); % ampacity greater than or equal to Idc
minrat = min(rating); % smallest ampacity greater than or equal to Idc
gauge = gauge(minrat); % gauge of wire with above rated ampacity
th = min(gauge)*.5067*0.00155; % thickness of wire in in^2

% Reference for voltage drop calc: http://www.adamselectric.coop/wp-content/uploads/2015/02/Voltage-Drop.pdf
drop = sqrt(3)*K*Idc*15/(gauge*1000); %Voltage Drop over 15 ft of cable

Ns = Ns + ceil(drop./DVdc); % Adds number of cells to account for voltage drop

%% Size up Battery Packs based on required energy 
E = DC.*DVdc.*Ns.*Np; % Battery Pack Energy in [kWh]

% Increases the number of the batteries in each pack to meet energy
% requiement by increasing number of batteries in parallel
Np = (E<Ebatt).*ceil(((Ebatt./E-1).*Np))+Np;

%% Mass of each Battery Pack

%Call tank function for tank length [in]
total_length = intertank.h + tank(1).dome_h*2+tank(1).h + tank(2).h;
cable = 0.1*th*total_length % mass of the 15 ft cable in the raceway
battmass  = 1.1*Ns.*Np.*Dmass; % Mass of each battery pack with a safety factor of 1.1

mass = cable+battmass; %mass of cable plus batteries(with a safety factor of 10% added to battery mass)
%% Select Battery Pack based on mass

[mass,index] = min(mass); % minimum battery pack mass
batt = struct('Name',Dname(index),'Ns',Ns(index),'Np',Np(index)); % Parameters of least massive battery pack

end