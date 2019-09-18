clc;
clear;
% First Iteration Setup
warning('off', 'all');
load_system('Purdue_Sim');
Pe = 5; %[psi]
L_star = 40;
CR = 5;
Y = 75000; % [for copper]
rho = 0.3190; %[for copper]
eta_cstar = 0.9;
eta_cf = 0.95;
prop = 1;
ffc = 0.1;
diameter = 18; %[in]
diameter_ft = diameter / 12; %[ft]
rocket_length = 24; %[ft]
i = 1;
j = 1;
altitude = [0];
ISP = 250;
load_system('Purdue_Sim');
lambda = 0.6;
mass_propellant = 200000/ISP;
mass_inert = mass_propellant*(1/lambda - 1);
mass_payload = 5;
wet_mass = mass_propellant + mass_inert + mass_payload; %[lbs]
Thrust = 5000; %[lbs]
Throttled_Thrust = Thrust;
Throttle_start_time = 0;
Throttle_end_time = 0;
eul_0 = [0,0,0]; %initial orientation [rad]
k_quat = 1; 
quat_statename = 'quat';
initial_height = 0; %[ft]
data = load('aerodata.csv');
m_dot = -Thrust/ISP; %[lbs/sec]
m_dot_slugs = m_dot / 32.2; %[slugs/sec]
Throttled_m_dot = m_dot_slugs;
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
max_q = max(q)/144;
%Optimization Loop
n = 10;
for prop = 1
    j = 1;
    if prop == 1
        D.fuel = 'methane';
        name = 'LNG';
        fuel_density = 24.70619;
    end
    if prop == 2
        D.fuel = 'rp-1';
        name = 'RP-1';
        fuel_density = 50.57;
    end
    if prop == 3
        D.fuel = 'ethanol';
        name = 'Ethanol';
        fuel_density = 49.25;
    end   
    for Thrust = 5000
        i = 1;
        for Pc = linspace(500, 850, n)
            k = 1;
            for OF = linspace(2, 3, n) 
                [mass_engine,ISP,isp_vac,cf,cf_vac,Dt,Dc,De,m_dot] = engineMass(OF, Pc, Pe, L_star, CR, Y, rho, eta_cstar, eta_cf, Thrust, prop, ffc); %engine performance function Ryan Strelau
                mass_propellant_0 = 200000/ISP; %Cal
                burn_time = 200000/ Thrust; %[sec]
                mass_propellant = mass_propellant_0 + mass_propellant_0*0.02 + mass_propellant_0*0.08/(1+OF) + mass_propellant_0*0.05/(1+1/OF); 
                [propellant, material] = set_tank_param(mass_propellant/(1+OF), mass_propellant/(1+1/OF), name, fuel_density);
                max_hydro_p = hydro_accel_calc(acceleration, oxidizer_m_array, diameter);        
                tank = tank_mass(propellant, 50, diameter, material(1), max_q, max_hydro_p);
                D.Pc = Pc;
                D.ox = 'oxygen';
                D.m_dot_ox = m_dot/(1+1/OF);
                D.m_dot_fuel = m_dot/(1+OF);
                D.tb = burn_time;
                mass_pump = PumpSizingModel(D);
                %mass_pump = 0;
                mass_sys = 100;
                mass_inert = mass_engine + mass_pump + tank(1).mass +tank(2).mass + mass_sys+50;
                tank_weight = tank(1).mass + tank(2).mass;
                mass_payload = 5;
                wet_mass = mass_propellant + mass_inert + mass_payload;%[lbs]
                dry_mass = wet_mass - m_dot * burn_time; %[lbs]
                [t, h] = ode45(@(t, h) vertical_launch(t, h, Thrust, wet_mass, m_dot), [0, 200], [0, 0]);
                dhdt = vertical_launch(t, h, Thrust, wet_mass, m_dot);
                height_max_OF(prop, j, i, k) = max(h(:, 1));
                P_c(prop, j, i, k) = Pc;
                m_dot_o(prop, j, i, k) = D.m_dot_ox;
                m_dot_f(prop, j, i, k) = D.m_dot_fuel;
                O_to_F(prop, j, i, k) = OF;
                engine(prop, j, i, k) = mass_engine;
                mass_fuel(prop, j, i, k) = mass_propellant/(1+OF);
                mass_ox(prop, j, i, k) = mass_propellant/(1+1/OF);
                m_tank(prop, j, i, k) = tank_weight;
                m_pumps(prop, j, i, k) = mass_pump;
                tank_volume_fuel(prop, j, i, k) = tank(1).volume;
                tank_volume_lox(prop, j, i, k) = tank(2).volume;
                tank_heights_fuel(prop, j, i, k)= tank(1).h;
                tank_heights_lox(prop, j, i, k) = tank(2).h;
                tank_thickness(prop, j, i, k) = tank.thick;
                tank_mat_v_fuel(prop, j, i, k) = tank(1).mat_v;
                tank_mat_v_lox(prop, j, i, k) = tank(2).mat_v;
                tank_mass_fuel(prop, j, i, k) = tank(1).mass;
                tank_mass_lox(prop, j, i, k) = tank(2).mass;
                ISP_output(prop, j, i, k) = ISP;
                wet_mass_ouput(prop, j, i, k) = wet_mass;
                dry_mass_output(prop, j, i, k) = dry_mass;
                k = k + 1;
            end
            height_max_Pc(prop, j, i) = max(height_max_OF(prop, j, i, :));
            i = i + 1;
        end
        height_max_Thrust(prop, j) = max(height_max_Pc(prop, j,:));
        j = j+1;
    end
    maximum = max(max(max(height_max_OF)));
    index = find(height_max_OF == maximum);
    wet_mass_ouput = reshape(wet_mass_ouput, [100, 1]);
    dry_mass_output = reshape(dry_mass_output, [100, 1]);
    tank_volume_fuel = reshape(tank_volume_fuel, [100, 1]);
    tank_volume_lox = reshape(tank_volume_lox, [100, 1]);
    tank_thickness = reshape(tank_thickness, [100, 1]);
    tank_mat_v_fuel = reshape(tank_mat_v_fuel, [100, 1]);
    tank_mat_v_lox = reshape(tank_mat_v_lox, [100, 1]);
    tank_mass_fuel = reshape(tank_mass_fuel, [100, 1]);
    tank_mass_lox = reshape(tank_mass_lox, [100, 1]);
    tank(1).volume = tank_volume_fuel(index);
    tank(2).volume = tank_volume_lox(index);
    tank(1).h = tank_heights_fuel(index);
    tank(2).h = tank_heights_lox(index);
    tank(1).thick = tank_thickness(index);
    tank(2).thick = tank_thickness(index);
    tank(1).mat_v = tank_mat_v_fuel(index);
    tank(2).mat_v = tank_mat_v_lox(index);
    tank(1).mass = tank_mass_fuel(index);
    tank(2).mass = tank_mass_lox(index);
    fprintf('Wet Mass = %.2f lbm\n', wet_mass_ouput(index));
    fprintf('Dry Mass = %.2f lbm\n', dry_mass_output(index));
    Pc_output = reshape(P_c, [100, 1]);
    Pc_output(index)
end




        
        