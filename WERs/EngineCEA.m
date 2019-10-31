function [cstar, isp, M, gamma, T, rho, mu, Pr, Mw, k] = EngineCEA(Pc, OF)
addpath('../WERs/cea');
% CEA_ROCKET_EXAMPLE: Uses MATLAB CEA wrapper. For in-depth
% documentation read the headers of cea_rocket_run.m,
% cea_rocket_run_single.m, and cea_rocket_read.m

% Change this variable to true to rerun CEA instead of using saved values
CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;
inp('type') = 'eq';                   % Sets the type of CEA calculation
inp('p') = Pc;                        % Chamber pressure
inp('p_unit') = 'psi';                % Chamber pressure units
inp('o/f') = OF;                      % Mixture ratio
% inp('sup') = 70;                    % Supersonic area ratios
inp('pip') = [1.75, Pc/12];           % Pressure ratios
inp('fuel') = 'CH4(L)';               % Fuel name from thermo.inp
inp('fuel_t') = 111.64;               % Fuel inlet temperature
inp('ox') = 'O2(L)';                  % Ox name from thermo.inp
inp('ox_t') = 90.17;                  % Ox inlet temperature
inp('file_name') = 'EngineCEA.inp';   % Input/output file name
if CEA_RUN
    data = cea_rocket_run(inp);       % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

% The output data structure, called 'data' in this case, is also a MATLAB
% map. 'data' contains a single entry for each of the CEA calculation types
% listed ('eq' and 'fr'). For instance, if only 'fr' is listed, then 'data'
% will only contain a single entry under data('fr').
data_eq = data('eq');

% Use keys(data_eq) or keys(data_fr) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.

cstar = squeeze(data_eq('cstar'));
M = squeeze(data_eq('mach'));
gamma = squeeze(data_eq('gammas'));
T = squeeze(data_eq('t'));
rho = squeeze(data_eq('rho'));
mu = squeeze(data_eq('visc'));
Pr = squeeze(data_eq('prandtl'));
Mw = squeeze(data_eq('m'));
k = squeeze(data_eq('k'));
isp = squeeze(data_eq('isp'));

cstar = cstar(1)*3.281; % Converts [m/s] to [ft/s]
T = T*1.8; % Converts [K] to [R]
rho = rho*3.613E-5; % Converts [kg/m^3] to [lbm/in^3]
mu = mu*1.450E-4; % Converts [Pa-s] to [psi-s]
Mw = Mw*2.205; % Converts [kg/kmol] to [lbm/kmol]
k = k*0.5782; % Converts [W/m-K] to [Btu/hr-ft-R]

end