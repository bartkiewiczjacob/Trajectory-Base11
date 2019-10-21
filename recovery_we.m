function mass = recovery_we(m_land, v_land)

% Constants and Assumptions
% Requirements
req.v_land = v_land;

% General
const.m_land = m_land; %[lb] landing weight (= dry mass)
const.landing_alt = 5000; %[ft]
const.air.rho0 = atmos(const.landing_alt, 'units', 'US'); %[slug/ft^3] air density at sea level

% Main Parachute
const.main.C_D0 = 0.9; % drag coefficient of main parachute (see Knacke fig. 5-25)
const.main.m_cloth = 0.0115; %[lb/ft^2] specific mass of main parachute cloth
const.main.susp_N = 20; % number of suspension lines
const.main.susp_m = 0.0075; %[lb/ft] specific mass of suspension lines
const.main.susp_F = 1000; %[lb] strength of suspension line
const.main.susp_LD = 1; % ratio of suspension line length to parachute diameter
const.main.rtapes_NG = 20; % Number of gores
const.main.rtapes_m = const.main.susp_m; %[lb/ft] specific weight of radial tapes (similar to suspension lines)
const.main.rtapes_F = const.main.susp_F; %[lb] strength of radial tapes
const.main.deploybagratio = 0.05; %[1] ratio of deployment bag mass on parachute mass

% Conversions
conv.lb2kg = 0.454; %[kg/lbm]
conv.kgm32slugft3 = 1/515.379;

%% Sizing

%%%%% Main Parachute
% Parachute size and weight
calc.main.q = 1/2 * const.air.rho0 * req.v_land^2; % Dynamic pressure at sea level
calc.main.CDS0 = const.m_land / calc.main.q; % Parachute drag area
calc.main.S0 = calc.main.CDS0/const.main.C_D0; %[ft^2] Parachute surface area
calc.main.D0 = sqrt(4*calc.main.S0 / pi); %[ft] Parachute nominal diameter
calc.main.A_cloth = 2*pi*(calc.main.D0/2)^2; %[ft^2] Surface area of parachute cloth assuming a hemisphere
calc.main.m_cloth = const.main.m_cloth * calc.main.A_cloth; %[lb] parachute cloth weight

% Suspension lines
calc.main.susp_L = calc.main.D0 * const.main.susp_LD; %[ft] suspension line length
calc.main.m_susp = const.main.susp_N * calc.main.susp_L * const.main.susp_m * const.main.susp_F/1000; %[lb] suspension line weight

% Radial tapes
calc.main.m_rtapes = calc.main.D0/2 * const.main.rtapes_NG * const.main.rtapes_m * const.main.rtapes_F/1000; %[lb] mass of radial tapes

% Deployment bag
calc.main.m_deploybag = (calc.main.m_cloth + calc.main.m_susp + calc.main.m_rtapes) * const.main.deploybagratio;

% Total main parachute weight
calc.main.m_total = calc.main.m_cloth + calc.main.m_susp + calc.main.m_rtapes + calc.main.m_deploybag; %[lb] total mass is the sum of all sub masses


%%%%% Drogue
calc.drogue.D0 = 8; %[ft]
calc.drogue.m_cloth = 5; %[lb]
calc.drogue.susp_L = calc.drogue.D0 * const.main.susp_LD; %[ft]
calc.drogue.m_susp = const.main.susp_N * calc.drogue.susp_L * const.main.susp_m * const.main.susp_F/1000; %[lb] suspension line weight
calc.drogue.m_riser = 1; %[lb]
calc.drogue.m_deploybag = (calc.drogue.m_cloth + calc.drogue.m_susp + calc.drogue.m_riser) * const.main.deploybagratio;
calc.drogue.m_total = calc.drogue.m_cloth + calc.drogue.m_susp + calc.drogue.m_riser + calc.drogue.m_deploybag; %[lb]


%%%%% Other weights
calc.ejection.m_total = 4; %[lb]
calc.heatshield.m_total = 1; %[lb]


%%%%% Addition
mass = calc.main.m_total + calc.drogue.m_total + calc.ejection.m_total + calc.heatshield.m_total; %[lb]

end