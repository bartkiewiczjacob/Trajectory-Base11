clc;
clear;
warning('off', 'all');
load_system('Purdue_Sim');
diameter = 18; %[in]
diameter_ft = diameter / 12; %[ft]
rocket_length = 24; %[ft]
mass_propellant = 1250 - 500;
mass_inert = 500;
mass_payload = 5;
wet_mass = mass_propellant + mass_inert + mass_payload; %[lbs]
Thrust = 5000; %[lbs]
throttle = 0;
Throttled_Thrust = Thrust*(1 - throttle);
Throttle_start_time = 15;
Throttle_end_time = 35;
Impulse_Savings = (Thrust - Throttled_Thrust)*(Throttle_end_time - Throttle_start_time);
eul_0 = [0,0,0]; %initial orientation [rad]
k_quat = 1; 
quat_statename = 'quat';
initial_height = 0; %[ft]
data = load('aerodata.csv');
burn_time = (200000+Impulse_Savings)/Thrust; %[sec]
m_dot = mass_propellant/burn_time; %[lbs/sec]
m_dot_slugs = m_dot / 32.2; %[slugs/sec]
Throttled_m_dot = m_dot_slugs*(1 - throttle);
rocket_length = 24; %[ft]
dry_mass = wet_mass - m_dot * burn_time; %[lbs]
wet_mass_slugs = wet_mass / 32.2; %[slugs]
dry_mass_slugs = dry_mass / 32.2; %[slugs]
I_wet = 1/12 * wet_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                   0, rocket_length^2 + 3*(diameter_ft/2)^2, 0; ...
                   0, 0, rocket_length^2 + 3*(diameter_ft/2)^2]; %[wet mass Inertia Tensor]
I_dry = 1/12 * dry_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                   0, rocket_length^2 + 3*(diameter_ft/2)^2, 0; ...
           0, 0, rocket_length^2 + 3*(diameter_ft/2)^2]; %[dry mass Inertia Tensor]
sim('Purdue_Sim'); % runs the simulation