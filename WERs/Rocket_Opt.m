%% ROCKET OPTIMIZATION CODE

clear
clc

% Initial Values
len_nose = 34.5; % in (plz change)
len_bottom = 12; % in (plz change)
m_nose = 6; % lbm (plz change)
m_rec = 2; %lbm (plz change)
m_plumbing = 15; % lbs (plz change)

%% INCREMENT VARIABLES

p_start = 500; % psi
p_inc = 50; % psi
p_end = 500; % psi
T_start = 1200; % lbf
T_end = 1200; % lbf
T_inc = 50; % lbf
OF = 3; % oxidizer to fuel ratio

% Altitude requirements
min_alt_goal = 42500; % ft
max_alt_goal = 47500; % ft

pressures = p_start:p_inc:p_end;
avail_inner_diameters = [5.25 5.5 5.75 6 6.5 7 7.5]; % in
thrusts = T_start:T_inc:T_end;

%% INITIALIZE RESULT MATRICES
% tank pressure = rows, inner diameter = columns, thrust = up
dry_mass = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
mass_fins = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
mass_eng = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
mass_tank = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
mass_str = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
mass_avionics = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
rocket_length = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
tank_length = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
str_length = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
Isp = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
outer_diameter = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
thickness_tank_CH4 = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
thickness_tank_LOX = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
thickness_fins = zeros(length(pressures), length(avail_inner_diameters), length(thrusts), 4);
max_alt = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
t_apo = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
max_vel = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
max_Mach = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
max_acc = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
alt_max_vel = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
time_vec = zeros(length(pressures), length(avail_inner_diameters), length(thrusts), 75000);
velocity_vec = zeros(length(pressures), length(avail_inner_diameters), length(thrusts), 75000);
altitude_vec = zeros(length(pressures), length(avail_inner_diameters), length(thrusts), 75000);
CH4_tank_outer = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
thrust2weight = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
heat_flux = zeros(length(pressures), length(avail_inner_diameters), length(thrusts));
prop_outputs = zeros(length(pressures), length(avail_inner_diameters), length(thrusts), 12);

%% ITERATION VARIABLES
row = 1;
col = 1;
up = 1;
i = 1;

anySuccess = 0;

%% CALCULATIONS
% tic
for thrust_eng = thrusts
    row = 1;
    for tank_pressure = pressures
        col = 1;
        for inner_diameter = 5.75
            
            index = [row col up]
            
            [m_tank, o_d, odCH4, tLOX, tCH4, len_tank] = tankWER(inner_diameter, tank_pressure);

            [m_eng,  Isp_eng, T2W, q_t, output ] = EngineSizing_FilmCooling(thrust_eng, tank_pressure, o_d, OF);

            [~, m_str, len_str] = Vehicle_WER(o_d, 0);

            len_tot = len_nose + len_str + len_tank + len_bottom;

            [m_fins, t_fins] = WER_Fins(len_tot, o_d);
            
            [m_avionics] = avionicsWER(2, 2, 30);

            % [m_rec] = recWER();

            m_tot = m_eng + m_fins + m_tank + m_rec + m_str + m_nose + m_avionics;
            
%             Isp_eng
%             m_tot*1.3
%             thrust_eng
%             o_d
                        
            [alt, t, v_max, v_max_alt, Mach_num_max, acc_max, velocity, altitude, time]...
                = Function_1DOF(o_d, thrust_eng, m_tot*1.3, Isp_eng);
            
%             alt
            
            if alt > min_alt_goal && alt < max_alt_goal
            success(i,:) = [row col up];
            i=i+1;
            anySuccess = anySuccess + 1;
            end

            % Result matrices of all values
            dry_mass(row, col, up) = m_tot;
            mass_fins(row, col, up) = m_fins;
            mass_eng(row, col, up) = m_eng;
            mass_tank(row, col, up) = m_tank;
            mass_str(row, col, up) = m_str;
            mass_avionics(row, col, up) = m_avionics;
            rocket_length(row, col, up) = len_tot;
            tank_length(row, col, up) = len_tank;
            str_length(row, col, up) = len_str;
            Isp(row, col, up) = Isp_eng;
            outer_diameter(row, col, up) = o_d;
            thickness_tank_CH4(row, col, up) = tCH4;
            thickness_tank_LOX(row, col, up) = tLOX;
            thickness_fins(row, col, up, :) = t_fins;
            max_alt(row, col, up) = alt;
            t_apo(row, col, up) = t;
            max_vel(row, col, up) = v_max;
            max_Mach(row, col, up) = Mach_num_max;
            max_acc(row, col, up) = acc_max;
            alt_max_vel(row, col, up) = v_max_alt;
            time_vec(row, col, up, :) = time;
            altitude_vec(row, col, up, :) = altitude;
            velocity_vec(row, col, up, :) = velocity;
            thrust2weight(row, col, up) = T2W;
            heat_flux(row, col, up) = q_t;
            prop_outputs(row, col, up, :) = output;

            col = col + 1;
        end
        row = row + 1;
    end
    up = up + 1;
end

fprintf('Number of Successes: %i\n', anySuccess);
% fprintf('Simulation time = %.3f s\n', runtime);

%% COMPILE SUCCESSFUL COMBINATIONS FOR EXPORT

if anySuccess > 0
    
    [numSuccesses,~] = size(success);
    
    general_values = zeros(numSuccesses, 27);
    success_prop_outputs = zeros(numSuccesses, 12);
    success_fin_dims = zeros(numSuccesses, 4);
    
    for k = 1:numSuccesses
        
        general_values(k,:) = [ success(k,1), success(k,2), success(k,3),...
                pressures(success(k,1)),...
                avail_inner_diameters(success(k,2)),...
                thrusts(success(k,3)),...
                max_alt(success(k,1),success(k,2),success(k,3)),...
                dry_mass(success(k,1),success(k,2),success(k,3)),...
                mass_fins(success(k,1),success(k,2),success(k,3)),...
                mass_eng(success(k,1),success(k,2),success(k,3)),...
                mass_tank(success(k,1),success(k,2),success(k,3)),...
                mass_str(success(k,1),success(k,2),success(k,3)),...
                mass_avionics(success(k,1),success(k,2),success(k,3)),...
                rocket_length(success(k,1),success(k,2),success(k,3)),...
                tank_length(success(k,1),success(k,2),success(k,3)),...
                str_length(success(k,1),success(k,2),success(k,3)),...
                Isp(success(k,1),success(k,2),success(k,3)),...
                outer_diameter(success(k,1),success(k,2),success(k,3)),...
                thickness_tank_CH4(success(k,1),success(k,2),success(k,3)),...
                thickness_tank_LOX(success(k,1),success(k,2),success(k,3)),...
                t_apo(success(k,1),success(k,2),success(k,3)),...
                max_vel(success(k,1),success(k,2),success(k,3)),...
                max_Mach(success(k,1),success(k,2),success(k,3)),...
                max_acc(success(k,1),success(k,2),success(k,3)),...
                alt_max_vel(success(k,1),success(k,2),success(k,3)),...
                thrust2weight(success(k,1),success(k,2),success(k,3)),...
                heat_flux(success(k,1),success(k,2),success(k,3)),...
           ];
               
        success_prop_outputs(k,:) = prop_outputs(success(k,1),success(k,2),success(k,3));
        success_fin_dims(k,:) = thickness_fins(success(k,1),success(k,2),success(k,3));
               
    end
end

%% EXPORT

csvwrite('RESULTS P 400-800 T 800-1200.csv', general_values);
csvwrite('PROP OUTPUTS P 400-800 T 800-1200.csv', success_prop_outputs);
csvwrite('FIN DIMS P 400-800 T 800-1200.csv', success_fin_dims);