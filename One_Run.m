clc;
clear;
warning('off', 'all');
load_system('Purdue_Sim');
Cd_data = load('cd_M.csv');
Cn_data = load('cl_M_a.csv');
Cp_data = load('cp_M_a.csv');
alpha_data = load('alpha_data.csv');
M_data = load('M_data.csv');
diameter = 18; %[in]
diameter_ft = diameter/12; %[ft]
rocket_length = 24; %[ft]
S = pi*diameter_ft^2/4;
throttle = 0;
latitude = 32.9861;
longitude = -106.9717;
margin = 0;
wet_mass = 1444.12 + margin;%[lbs]
mass_inert = 523.36 + margin;
Thrust = 5100; %[lbs]
Throttled_Thrust = Thrust*(1 - throttle);
Throttle_start_time = 15;
Throttle_end_time = 35;
De = 8.5883;
Dt = 3.1665;
At = pi*Dt^2/4;
Pc = 450;
Pe = 10;
AR = De^2/Dt^2;
Impulse_Savings = (Thrust - Throttled_Thrust)*(Throttle_end_time - Throttle_start_time);
eul_0 = [0,0,0]; %initial orientation [rad]
k_quat = 1; 
quat_statename = 'quat';
initial_height = 4595; %[ft]
burn_time = (200000+Impulse_Savings)/Thrust; %[sec]
m_dot = 23.35; %[lbs/sec]
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
for i = 1
    cg_offset_r = normrnd(0, 0.1)/12;
    cg_offset_theta = rand()*2*pi;
    y_offset = cg_offset_r*sin(cg_offset_theta);
    z_offset = cg_offset_r*cos(cg_offset_theta);
    sim('Purdue_Sim'); % runs the simulation
    %max_q = max(q)/144;
    result(i) = max(ans.altitude)*0.0003048;
    i
end
histogram(result)
    