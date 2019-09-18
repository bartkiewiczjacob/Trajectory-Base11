clc;
clear;
warning('off', 'all');
load_system('Purdue_Sim');
diameter = 18; %[in]
diameter_ft = diameter / 12; %[ft]
rocket_length = 24; %[ft]
i = 1;
j = 1;
altitude = 999999999999;
margin = 160;
for throttle = 0
    margin = 160;
    altitude = 999999999999;
    while max(altitude) > 410105
        mass_propellant = 1250.19 - 495.81;
        mass_inert = 495.81 - 100 + margin;
        mass_payload = 5;
        wet_mass = mass_propellant + mass_inert + mass_payload; %[lbs]
        Thrust = 5000; %[lbs]
        Throttled_Thrust = Thrust*(1 - throttle);
        Throttle_start_time = 15;
        Throttle_end_time = 35;
        Impulse_Savings = (Thrust - Throttled_Thrust)*(Throttle_end_time - Throttle_start_time);
        eul_0 = [0,0,0]; %initial orientation [rad]
        k_quat = 1; 
        quat_statename = 'quat';
        initial_height = 4595; %[ft]
        data = load('aerodata.csv');
        burn_time = (200000+Impulse_Savings)/Thrust; %[sec]
        m_dot = 14.14+4.71; %[lbs/sec]
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
        max_q = max(q)/144;
        margin = margin + 1;
    end
%     margin_output(i) = margin;
%     i = i + 1;
%     figure(1)
%     plot(time, q);
%     hold on;
%     title('q vs. time');
%     ylabel('q [psf]');
%     xlabel('time [sec]');
end
% throttle = linspace(0, 0.3, 10);
% legendCell = cellstr(num2str(throttle', 'Throttle Percentage =%-d'));
% legend(legendCell)
% figure(2)
% plot(throttle, margin_output);
% title('Non-Fuel System Mass vs. Throttle Percentage');
% xlabel('Throttle Percentage');
% ylabel('Inert Mass [lbm]');
maximum = max(q);
index = find(q==maximum);
density(index)
velocity(index)