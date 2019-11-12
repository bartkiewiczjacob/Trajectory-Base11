clc;
clear;
% First Iteration Setup
warning('off', 'all'); %Turns off Simulink Warnings to clear up command.
load_system('Purdue_Sim'); %Preloads Sim to decrease run time
mat_prop = readtable('Material_prop.xlsx'); %Loads material properties for airframe
Pe = 10; %Exit pressure for nozzle[psi]
De = 8.5883;
Dt = 3.1665;
At = pi*Dt^2/4;
Pc = 450;
AR = De^2/Dt^2;
L_star = 40; %Engine parameter
CR = 5; %Engine Parameter
Y = 75000; % Engine parameter [for copper]
rho = 0.3190; %Engine parameter[for copper]
eta_cstar = 0.9; 
eta_cf = 0.95;
ffc = 0.1;
diameter = 18; %Rocket Diameter[in]
diameter_ft = diameter / 12; %Rocket Diameter[ft]
rocket_length = 24; %Rocket Length [ft]
% Initial Run Parameters
ISP = 190.48; %ISP [sec]
lambda = 0.6; %Propellant Mass Fraction
OF = 1.45; %Oxygen to Fuel Ratio
mass_propellant = 200000/ISP; %Mass of Propellant [lbm]
mass_inert = mass_propellant*(1/lambda - 1); %Inert Mass [lbm]
mass_payload = 5; %Payload Mass [lbm]
wet_mass = mass_propellant + mass_inert + mass_payload; %GLOW [lbm]
Thrust = 5100; %[lbs]
% Throttle Info For 6DOF (Don't Change This)
Throttled_Thrust = Thrust;
Throttle_start_time = 0;
Throttle_end_time = 0;
%Further Simulink Initialization
eul_0 = [0,0,0]; %initial orientation [rad]
k_quat = 1; 
quat_statename = 'quat';
initial_height = 4595; %[ft]
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
latitude = 32.9861;
longitude = -106.9717;
S = pi*diameter_ft^2/4;
Cd_data = load('cd_M.csv');
Cn_data = load('cl_M_a.csv');
Cp_data = load('cp_M_a.csv');
alpha_data = load('alpha_data.csv');
M_data = load('M_data.csv');
sim('Purdue_Sim'); % runs the simulation
max_q = max(q)/144; %Initial Max Q [psi]
Cd = 0.3; %Drag Coefficient
S = pi*(0.0254*18)^2/4; %Reference Area [m^2]
Aerodata = load('mach_cd.csv');
%Optimization Loop
i = 1;
j = 1;
n = 10; %number of points in optimization 
pressfed = 0; % 1 if pressure fed 0 if pump fed
tank_temp = -150; %Tank temperature at empty
x = 1;
margin = 0;
for diameter = 16
    S = pi*(0.0254*diameter)^2/4;
    for prop = 3 %Propellant Choice Loop
        j = 1; %Follow existing format for propellant data (Pumps will need to be updated sepereately to accept new fuels)
        if prop == 1
            D.fuel = 'methane';
            name = 'LNG';
            fuel_density = 24.70619;
            Pc_floor = 750;
            OF_floor = 2;
        end
        if prop == 2
            D.fuel = 'rp-1';
            name = 'RP-1';
            fuel_density = 50.57;
            Pc_floor = 350; %Change to 200 for Pressure Fed. Otherwise 300
            OF_floor = 2;
        end
        if prop == 3
            D.fuel = 'ethanol';
            name = 'Ethanol';
            fuel_density = 49.25;
            Pc_floor = 350;
            OF_floor = 1;
        end
        for Thrust = 5100 %Thrust Loop 7000 max for pressure fed, 6000 max for pump
            i = 1;
            for Pc = 450 %Chamber Pressure Loop
                k = 1;
                for OF = 1.45  %Oxygen to Fuel Ratio Loop
                %Rocket Design
                    [mass_engine,ISP,isp_vac,cf,cf_vac,Dt,Dc,De,m_dot] = engineMass(OF, Pc, Pe, L_star, CR, rho, eta_cstar, eta_cf, Thrust, prop, ffc); %engine performance function Ryan Strelau & Matt Schmale
                    OF = OF/1.1;
                    mass_propellant_0 = 200000/ISP; %Propellant Mass [lbm]
                    burn_time = 200000/ Thrust; %[sec]
                    mass_propellant = mass_propellant_0 + mass_propellant_0*0.02 + mass_propellant_0*0.08/(1+OF) + mass_propellant_0*0.05/(1+1/OF); %Updates Propellant Mass for pump leakage and ullage
                    [tank, helium, material,max_hydro_p] = calc_tanks(Pc,pressfed,tank_temp,diameter,max_q, mass_propellant/(1+OF), mass_propellant/(1+1/OF), oxidizer_m_array, acceleration, fuel_density); %Calculates tank weights and dimensions Ben Davis
                    intertank = intertank_mass(tank, material(1)); %Calculates weight and dimensions of structure between tanks Ben Davis
                    L1 = 2*diameter; %Engine and Pump Structure Length
                    fin = fin_analysis(diameter, 3, mat_prop); %Calculates Fin dimensions and mass Shushrut Prahbu
                    boattail = Boattail(diameter, L1, De,fin(2)); %Calculates Boatail 
                    %Pump Structure Definition
                    D.rpm = 45000;
                    D.Pc = Pc; %Chamber Pressure (defined from loop)
                    D.ox = 'oxygen'; %Oxidizer Name (don't change)
                    D.m_dot_ox = m_dot/(1+1/OF); %Mass flow of Oxygen [lbm]
                    D.m_dot_fuel = m_dot/(1+OF); %Mass flow of Fuel [lbm]
                    D.tb = burn_time; %Burn Time [sec]
                    mass_pump = PumpSizingModel(D); %Determines mass of pumps Kyle Runkle 
                    %mass_pump = 0; %Uncomment if Pressure Fed
                    mass_avionics = 30; %Avionics Mass Humza Nasir
                    mass_RCS = 20; %RCS Mass
                    mass_inert = mass_engine + mass_pump + tank(1).mass +tank(2).mass + intertank.mass + mass_avionics + mass_RCS + boattail(1)+4*fin(1)+helium.total_mass + 10+margin;
                    tank_weight = tank(1).mass + tank(2).mass + intertank.mass;
                    mass_payload = 5;
                    wet_mass = mass_propellant + mass_inert + mass_payload;%[lbs]
                    dry_mass = wet_mass - m_dot * burn_time; %[lbs]

                    mass_recovery = recovery_we(dry_mass, 30);
                    dry_mass = (dry_mass + mass_recovery);
                    mass_margin = dry_mass*0.3;
                    wet_mass = wet_mass + mass_recovery + mass_margin;
                    dry_mass = dry_mass + mass_margin;
                    lambda = (wet_mass - dry_mass)/wet_mass;
                %Simulation of Flight
                    [t, h] = ode45(@(t, h) vertical_launch(t, h, Thrust, m_dot, burn_time, Cd, De, Dt, Pc, Pe, S), [0 300], [1400, 0, wet_mass*0.453592]);
                    max_h = max(h(:, 1));
                    index_h = find(h(:, 1) == max_h);
                    h = h(1:index_h, :);
                    t = t(1:index_h, :);
                    [temp, A, P, RHO] = atmoscoesa(h(:,1));
                    G = h(:,3)*9.8; %[N]
                    q = 0.5*RHO.*h(:,2).*h(:,2);
                    Mach = h(:,2)./A;
                    max_q = max(q)*0.000145038;
                    index_q = find(q == max_q/0.000145038);
                    h(index_q, 1);
                    Thrust_adj = getEngineAtAlt(Thrust, De^2/Dt^2, P*0.000145038, Pc, pi*Dt^2/4, Pe);
                    thrust_acceleration = Thrust_adj*4.45./h(:,3);
                    drag_acceleration = -Cd.*q*S./h(:,3);
                    acceleration = (1./h(:,3)).*(-Cd.*q*S + Thrust_adj*4.45.*(t < 40) - G);
                    %acceleration = acceleration*3.28/32.2;
                    oxidizer_m_array = h(:,3)/(1+1/OF)*2.20462;
                %Design Parameter Recording
                    De_output(x, prop, j, i, k) = De;
                    Dt_output(x, prop, j, i, k) = Dt;
                    tb(x, prop, j, i, k) = burn_time;
                    P_c(x,prop, j, i, k) = Pc;
                    m_dot_o(x, prop, j, i, k) = D.m_dot_ox;
                    m_dot_f(x, prop, j, i, k) = D.m_dot_fuel;
                    m_dot_output(x, prop, j, i, k) = m_dot;
                    O_to_F(x, prop, j, i, k) = OF;
                    engine(x, prop, j, i, k) = mass_engine;
                    mass_fuel(x, prop, j, i, k) = mass_propellant/(1+OF);
                    mass_ox(x, prop, j, i, k) = mass_propellant/(1+1/OF);
                    m_tank(x, prop, j, i, k) = tank_weight;
                    m_pumps(x, prop, j, i, k) = mass_pump;
                    tank_volume_fuel(x, prop, j, i, k) = tank(1).volume;
                    tank_volume_lox(x, prop, j, i, k) = tank(2).volume;
                    tank_heights_fuel(x, prop, j, i, k)= tank(1).h;
                    tank_heights_lox(x, prop, j, i, k) = tank(2).h;
                    tank_thickness(x, prop, j, i, k) = tank.thick;
                    tank_mat_v_fuel(x, prop, j, i, k) = tank(1).mat_v;
                    tank_mat_v_lox(x, prop, j, i, k) = tank(2).mat_v;
                    tank_mass_fuel(x, prop, j, i, k) = tank(1).mass;
                    tank_mass_lox(x, prop, j, i, k) = tank(2).mass; 
                    ISP_output(x, prop, j, i, k) = ISP;
                    wet_mass_output(x, prop, j, i, k) = wet_mass;
                    dry_mass_output(x, prop, j, i, k) = dry_mass;
                    m_recovery(x, prop, j, i, k) = mass_recovery;
                    thrust_output(x, prop, j, i, k) = Thrust;
                    %drag_loss(prop, j, i, k) = trapz(t, drag_acceleration);
                    gravity_loss(x, prop, j, i, k) = -32.2*burn_time;
                    helium_mass(x, prop, j, i, k) = helium.total_mass;
                    mass_prop(x, prop, j, i, k) = mass_propellant;
                    maxQ(x, prop, j, i, k) = max_q;
                    maxMach(x, prop, j, i, k) = max(Mach);
                    tfinal(x, prop, j, i, k) = t(end);
                    lambda_output(x, prop, j, i, k) = lambda;
                    if Thrust/(wet_mass) < 0
                        height_max_OF(x, prop, j, i, k) = 0;
                    else
                        height_max_OF(x, prop, j, i, k) = max(h(:, 1));
                    end
                    k = k + 1;
                end
                height_max_Pc(x, prop, j, i) = max(height_max_OF(x, prop, j, i, :));
                i = i + 1;
            end
            height_max_Thrust(x, prop, j) = max(height_max_Pc(x, prop, j,:));
            j = j+1;
        end
        height_max_Fuel(x, prop) = max(height_max_Thrust(x, prop, :));
        maximum = max(max(max(max(max(height_max_OF(x, prop, :, :, :))))));
        index(prop) = find(height_max_OF == maximum);
        tank(1).volume = tank_volume_fuel(index(prop));
        tank(2).volume = tank_volume_lox(index(prop));
        tank(1).h = tank_heights_fuel(index(prop));
        tank(2).h = tank_heights_lox(index(prop));
        tank(1).thick = tank_thickness(index(prop));
        tank(2).thick = tank_thickness(index(prop));
        tank(1).mat_v = tank_mat_v_fuel(index(prop));
        tank(2).mat_v = tank_mat_v_lox(index(prop));
        tank(1).mass = tank_mass_fuel(index(prop));
        tank(2).mass = tank_mass_lox(index(prop));
        fprintf('Results for %s\n', name);
        fprintf('Wet Mass = %.2f lbm\n', wet_mass_output(index(prop)));
        fprintf('Dry Mass = %.2f lbm\n', dry_mass_output(index(prop)));
        fprintf('Engine Mass = %.2f lbm\n', engine(index(prop)));
        fprintf('Tank Mass = %.2f lbm\n', m_tank(index(prop)));
        fprintf('Pump Mass = %.2f lbm\n', m_pumps(index(prop)));
        fprintf('Recovery Mass = %.2f lbm\n', m_recovery(index(prop)));
        fprintf('Chamber Pressure = %.2f psi\n', P_c(index(prop)));
        fprintf('ISP = %.2f s\n', ISP_output(index(prop)));
        fprintf('OF = %.2f\n', O_to_F(index(prop)));
        fprintf('Mass Flow = %.2f lbm\n', m_dot_o(index(prop)) + m_dot_f(index(prop)));
        fprintf('Thrust = %.2f lbf\n', thrust_output(index(prop)));
        fprintf('Apogee = %.2f km\n', maximum/1000);
        fprintf('Max Mach = %.2f\n', maxMach(index(prop)));
        fprintf('Max Q = %.2f psi\n', maxQ(index(prop)));
        fprintf('Time to Apogee = %.2f s\n', tfinal(index(prop)));
%         Thrust = linspace(3000, 7000, n);
%         figure(1)
%         plot(Thrust, height_max_Thrust(prop, :, :, :));
%         title('Thrust vs. Max Altitude for Pump Fed');
%         hold on; 
    end
    height_max_Diameter(x) = max(height_max_Fuel(x, :)); 
%     figure(1)
%     ylabel('Max Height [km]');
%     xlabel('Thrust [lbm]');
%     hold off;
    x = x+1;
end
% for i = 1:10
%     Pc_plot(i) = max(max(max(max(max(height_max_OF(:,:, :, i, :))))));
%     index(i) = find(height_max_OF == Pc_plot(i));
%     lambda_plot(i) = lambda_output(index(i));
% end
% Pc = linspace(Pc_floor, 800, 500);
% plot(Pc, lambda_plot)
OF = linspace(1.1, 1.6, 100);
index = find(height_max_OF ~= 0);

plot(OF, height_max_OF(index))


        
        