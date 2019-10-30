function [ M,  Isp_corrected, T2W, q_t, output ] = EngineSizing_FilmCooling( F, Pc, D, OF )
% Written by Robert Groome - Oct 26, 2019
% Film cooling added by Matt Hoeper - Oct 27, 2019
% Heat flux calcs to come later
% Estimates the mass of an engine. Also provides information for specific
% impulse and thrust to weight ratio
% Inputs:
% F = Engine Thrust [lbf]
% Pc = Chamber Pressure [psi]
% D = Diameter of Vehicle [in]
% OF = Oxidizer to Fuel Ratio of Engine
% Outputs:
% M = Mass of the Engine [lbm]
% Isp = Specific Impulse of Engine [s]
% T2W = Thrust to Weight Ratio
% q_t = Heat flux in the throat [Btu/hr-ft^2]
% output = a vector of other useful engine variables

% CONTROLS NOTES: change input from Cp to tank_p and move associated calc
% up in function for use in program, 

delete *.inp; % fixed CEA bug

%% Combustion Properties
[cstar, ~, M, gamma, T, rho, mu, Pr, Mw, k] = EngineCEA(Pc, OF); 

%% Film Combustion Properties
[cstar_film, isp_film, M_film, gamma_film, T_film, rho_film, mu_film, Pr_film, Mw_film, k_film] = EngineCEA(Pc, 0.1); 
Pe = 12; %Exit pressure, slightly less than atmospheric pressure at sea level
Pa = 14.7; %Atmospheric pressure at sea level
Itot = 9208; % Total impulse for a O-class engine [lbf-s]

Me = M(3);          % Exit Mach number
gamma_c = gamma(1); % Chamber gamma
gamma_t = gamma(2); % Throat gamma
gamma_n = gamma(3); % Nozzle gamma
T_c = T(1);         % Chamber temperature [R]
T_t = T(2);         % Throat temperature [R]
T_n = T(3);         % Nozzle temperature [R]
rho_c = rho(1);     % Chamber density [lbm/in^3]
rho_t = rho(2);     % Throat density [lbm/in^3]
rho_n = rho(3);     % Nozzle density [lbm/in^3]
mu_c = mu(1);       % Chamber viscosity [psi-s]
mu_t = mu(2);       % Throat viscosity [psi-s]
mu_n = mu(3);       % Nozzle viscosity [psi-s]
Pr_c = Pr(1);       % Chamber Prandtl number
Pr_t = Pr(2);       % Throat Prandtl number
R_c = 10.73/Mw(1);  % Chamber gas constant [ft^3-psi/lbm-R]
R_t = 10.73/Mw(2);  % Throat gas constant [ft^3-psi/lbm-R]
R_n = 10.73/Mw(3);  % Nozzle gas constant [ft^3-psi/lbm-R]
k_c = k(1);         % Chamber thermal conductivity [Btu/hr-ft-R]
k_t = k(2);         % Throat thermal conductivity [Btu/hr-ft-R]

%% Performance Evaluation
epsilon = 1/Me * ((2+(gamma_n-1)*Me^2)/(gamma_n+1))^((gamma_n+1)/(2*(gamma_n-1))); % Expansion ratio
cf = sqrt(((2*gamma_n^2)/(gamma_n-1)) * (2/(gamma_n+1))^((gamma_n+1)/(gamma_n-1)) * (1-(Pe/Pc)^((gamma_n-1)/gamma_n))); % Dimensionless jet thrust
cf = cf + (Pe/Pc - Pa/Pc)*epsilon; % Add dimensionless pressure thrust to get thrust coefficient
cf = cf*0.9; % 90% cf efficiency
cstar = cstar*0.80; % 90% c* efficiency [ft/s]
Isp = cf*cstar / 32.2; % Specific impulse [s]

%% Fluids Evaluation
De = D / 1.25; % Exit diameter assumed to be half of vehicle diameter [in]
At = pi*(De/2)^2 / epsilon; % Area of the throat [in^2]
Dt = 2*sqrt(At/pi); % Throat diameter [in]
mdot = F / Isp; % Required mass flow [lbf/s]
tb = Itot / F; % Burn time [s]
% mprop = mdot * tb; % Total propellant mass [lbf]
tank_p = Pc/(0.75*0.8);

%% Film Cooling
mdot_core_fuel = mdot/(1 + OF);
mdot_ox = mdot - mdot_core_fuel;
mdot_ffc = mdot_core_fuel * 0.15;
Isp_corrected = (mdot * Isp + mdot_ffc*isp_film(3)*0.8/32.2) / (mdot + mdot_ffc);

mdot_fuel_total = mdot_ffc + mdot_core_fuel;
mfuel = mdot_fuel_total * tb;
mox = mdot_ox * tb;
%% Engine Geometry and Mass Calc
CR = 5; % Contraction ratio
Ac = At * CR; % Chamber area [in^2]
Dc = 2*sqrt(Ac/pi); % Chamber diameter [in]
Lstar = 40; % Characteristic length (needs verification) [in]
Lc = Lstar / CR; % Length of chamber [in]
Ln = (D-Dt/2) / tan(deg2rad(15)); % Length of nozzle [in]
tc = 1.5 * (Pc*Dc/2) / (124000*(1.065-0.000222*T_c)); % Thickness of chamber wall using copper [in]
tn = 1.5 * (Pc*D/8) / (124000*(1.065-0.000222*T_n)); % Thickness of nozzle wall using copper [in]

V = pi*Dc*Lc*tc + pi*0.5*D*sqrt(Ln^2+(D/2)^2)*tn; % Volume of the chamber and nozzle [in^3]
M = V * 0.3; % Mass of the chamber and nozzle [lbm]
M = M + 5; % Estimating 5 lbs for injector and manifold [lbm]

T2W = F/M; % Thrust to weight ratio

%% Heat Flux Evaluation
Tr_c = T_c*(1+(gamma_c-1)/2*0.2^2*Pr_c^1/3); % Recovery temperature in the chamber [R]
Tam_c = 0.5*(T_c+Tr_c); % Arithmatic mean temperature in the chamber [R]
rhoam_c = interp1([T_n,T_t,T_c],[rho_n,rho_t,rho_c],Tam_c,'linear','extrap'); % Arithmatic mean density in the chamber [lbm/in^3]
muam_c = interp1([T_n,T_t,T_c],[mu_n,mu_t,mu_c],Tam_c,'linear','extrap'); % Arithmatic mean viscosity in the chamber [psi-s]

Tr_t = T_t*(1+(gamma_t-1)/2*Pr_t^1/3); % Recovery temperature at the throat [R]
Tam_t = 0.5*(T_t+Tr_t); % Arithmatic mean temperature in the throat [R]
rhoam_t = interp1([T_n,T_t,T_c],[rho_n,rho_t,rho_c],Tam_t,'linear','extrap'); % Arithmatic mean density in the chamber [lbm/in^3]
muam_t = interp1([T_n,T_t,T_c],[mu_n,mu_t,mu_c],Tam_t,'linear','extrap'); % Arithmatic mean viscosity in the chamber [psi-s]

Re_c = rho_c*sqrt(gamma_c*R_c*T_c/4636.8)*Dc / (32.2*mu_c); % Reynolds Number in the chamber
Re_t = rho_t*sqrt(gamma_t*R_t*T_t/4636.8)*Dt / (32.2*mu_t); % Reynolds Number in the throat
Nu_c = 0.026*Re_c^0.8*Pr_c^0.8*(rhoam_c/rho_c)^0.8*(muam_c/mu_c)^0.2; % Nustle Number in the chamber using Bartz Equation
Nu_t = 0.026*Re_t^0.8*Pr_t^0.8*(rhoam_t/rho_t)^0.8*(muam_t/mu_t)^0.2; % Nustle Number in the thorat using Bartz Equation
q_c = Nu_c*k_c/(12*Dc) * (Tr_c-T_c); % Heat flux in the chamber [Btu/hr-ft^2]
q_t = Nu_t*k_t/(12*Dt) * (Tr_t-T_t); % Heat flux in the throat [Btu/hr-ft^2]

%% Variable Collection
output = [mdot_ox, mdot_fuel_total, mfuel, mox, tank_p, Dt, epsilon, Lc, Ln, tc, tn, q_c]; % Collection of engine design parameters

end