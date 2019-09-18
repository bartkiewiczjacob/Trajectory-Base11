clc;
clear;
i = 1;
warning('on', 'all');
altitude = [0];
delta = 0.01;
for ISP = linspace(100, 380, 10)
    lambda_o = 0.1;
    altitude = [0];
    fprintf('Optimizing for ISP = %.2f s \n', ISP);
    while max(altitude) < 410105
        load_system('Purdue_Sim');
        lambda = lambda_o;
        mass_propellant = 200000/ISP;
        mass_inert = mass_propellant*(1/lambda - 1);
        mass_payload = 5;
        wet_mass = mass_propellant + mass_inert + mass_payload; %[lbs]
        Thrust = 5000; %[lbs]
        eul_0 = [0,0,0]; %initial orientation [rad]
        k_quat = 1; 
        quat_statename = 'quat';
        initial_height = 0; %[ft]
        data = load('aerodata.csv');
        m_dot = -Thrust/ISP; %[lbs/sec]
        m_dot_slugs = m_dot / 32.2; %[slugs/sec]
        diameter = 18; %[in]
        diameter_ft = diameter / 12; %[ft]
        rocket_length = 24; %[ft]
        burn_time = 200000/Thrust; %[sec]
        dry_mass = wet_mass + m_dot * burn_time; %[lbs]
        wet_mass_slugs = wet_mass / 32.2; %[slugs]
        dry_mass_slugs = dry_mass / 32.2; %[slugs]
        I_wet = 1/12 * wet_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, rocket_length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, rocket_length^2 + 3*(diameter_ft/2)^2]; %[wet mass Inertia Tensor]
        I_dry = 1/12 * dry_mass_slugs * [3*(diameter_ft/2)^2, 0, 0; ...
                           0, rocket_length^2 + 3*(diameter_ft/2)^2, 0; ...
                           0, 0, rocket_length^2 + 3*(diameter_ft/2)^2]; %[dry mass Inertia Tensor]
        sim('Purdue_Sim'); % runs the simulation
        lambda_o = lambda_o + delta;
        fprintf('Max Alt = %.0f ft\n', max(altitude));
    end
    min_lambda(i) = lambda_o - delta;
    i = i + 1;
end
sigma = 0.032851;
num_sigma = 4;
figure(1)
ISP = linspace(100, 380, 10);
x = 10^2:10^5;
y = 0.0158*log(x) + 0.7134;
y2 = y - sigma*num_sigma;
semilogx(x, y);
hold on;
semilogx(200000./ISP, min_lambda);
semilogx(x, y2);
hold off;
% plot(time, acceleration); 
title('Propellant Mass Fraction vs. Propellant Mass');
xlabel('Propellant Mass [lbm]');
ylabel('Propellant Mass Fraction');
grid on;






