
clc;
clear;
i = 1;
load_system('Purdue_Sim');
Thurst = 6500;
data = load('aerodata.csv');
m_dot = -22.7; %[lbs/sec]
m_dot_slugs = m_dot / 32.2; %[slugs/sec]
diameter = 18; %[in]
diameter_ft = diameter / 12; %[ft]
length = 28; %[ft]
burn_time = 200000/Thurst; %[sec]
prop_mass = -m_dot * burn_time;
wet_mass = prop_mass * 1.6; %[lbs]
dry_mass = wet_mass + m_dot * burn_time; %[lbs]
wet_mass_slugs = wet_mass / 32.2; %[slugs]
dry_mass_slugs = dry_mass / 32.2; %[slugs]
I_wet = 1/12 * wet_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, length^2 + 3*(diameter_ft/2)^2]; %[wet mass Inertia Tensor]
I_dry = 1/12 * dry_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, length^2 + 3*(diameter_ft/2)^2]; %[dry mass Inertia Tensor]
sim('Purdue_Sim');
disp(max(height));
disp(max(velocity));



