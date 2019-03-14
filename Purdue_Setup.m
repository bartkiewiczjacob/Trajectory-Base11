clc;
clear;
i = 1;
load_system('Purdue_Sim');
Thrust = 5000;
eul_0 = [0,0,0];
k_quat = 1;
quat_statename = 'quat';
initial_height = 4600;
data = load('aerodata.csv');
m_dot = -20; %[lbs/sec]
m_dot_slugs = m_dot / 32.2; %[slugs/sec]
diameter = 18; %[in]
diameter_ft = diameter / 12; %[ft]
rocket_length = 24; %[ft]
burn_time = 200000/Thrust; %[sec]
wet_mass = 1200; %[lbs]
dry_mass = wet_mass + m_dot * burn_time; %[lbs]
wet_mass_slugs = wet_mass / 32.2; %[slugs]
dry_mass_slugs = dry_mass / 32.2; %[slugs]
I_wet = 1/12 * wet_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, rocket_length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, rocket_length^2 + 3*(diameter_ft/2)^2]; %[wet mass Inertia Tensor]
I_dry = 1/12 * dry_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, rocket_length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, rocket_length^2 + 3*(diameter_ft/2)^2]; %[dry mass Inertia Tensor]
sim('Purdue_Sim');
j = 1;
Mach_op = [];
density_op = [];
for i = 1:1373
    for k = 1:length(Mach_op)
        if abs(Mach(i) - Mach_op(k)) < 0.001
            op(j) = time(i);
            j = j + 1;
        end
    end
    for x = 1:length(density_op)
        if abs(density(i) - density_op(x)) < 0.001
            op(j) = time(i);
            j = j + 1;
        end
    end
end


