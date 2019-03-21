clc;
clear;
i = 1;
load_system('Purdue_Sim');
Thrust = 5000;
eul_0 = [0,0,0];
k_quat = 1;
quat_statename = 'quat';
initial_height = 0;
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
Mach_op = [0.3 0.8 1.2 3.4];
density_op = [0.002377 0.002377/5 0.002377/10 0.002377/50 0.002377/100];
altitude_op = [];
n = 1;
m = 1;
for i = 1:1357
    for k = 1:length(Mach_op)
        if abs(Mach(i) - Mach_op(k)) < 0.0022
            mach_time_op(j) = time(i);
            j = j + 1;
        end
    end
    for x = 1:length(density_op)
        if abs(density(i) - density_op(x)) < 0.0000001
            density_time_op(n) = time(i);
            n = n + 1;
        end
    end
    for y = 1:length(altitude_op)
        if abs(altitude(i) - altitude_op(y)) < 10
            op(m) = time(i);
            m = m + 1;
        end
    end
end


