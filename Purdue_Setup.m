
clc;
clear;
i = 1;
load_system('Purdue_Sim');
Thrust = 5000;
initial_height = 4600;
data = load('aerodata.csv');
m_dot = -20; %[lbs/sec]
m_dot_slugs = m_dot / 32.2; %[slugs/sec]
diameter = 18; %[in]
diameter_ft = diameter / 12; %[ft]
length = 20; %[ft]
burn_time = 200000/Thrust; %[sec]
wet_mass = 1400; %[lbs]
dry_mass = wet_mass + m_dot * burn_time; %[lbs]
wet_mass_slugs = wet_mass / 32.2; %[slugs]
dry_mass_slugs = dry_mass / 32.2; %[slugs]
I_wet = 1/12 * wet_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, length^2 + 3*(diameter_ft/2)^2]; %[wet mass Inertia Tensor]
I_dry = 1/12 * dry_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, length^2 + 3*(diameter_ft/2)^2]; %[dry mass Inertia Tensor]
sim('Purdue_Sim')
%plot3(trajectory(:,1), trajectory(:,2), trajectory(:,3));
%grid on;
%view(30,30);
burn = trajectory - drift;


