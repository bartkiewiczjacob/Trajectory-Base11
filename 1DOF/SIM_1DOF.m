% Andrew Meyer
% 10/04/2019

% Script Initialization
clear, clc
% c_order = get(gca, 'ColorOrder');
% figNum = 0;


% Define constants
g_0 = 9.8; % m/s^2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Rocket Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dry mass
m_dry_lbm = 75; % lbm
m_dry = m_dry_lbm / 2.20462; % kg

% Drag coefficient
Cd = 0.4;

% Reference area
S = pi * (3*0.0254)^2; % m^2

% Thrust
thrust_lbf = 1000; % lbf
thrust = thrust_lbf * 4.44822; % N

% Specific Impulse
Isp = 185; % s

% Total Impulse (Max allowable = 9208 lbf*s)
impulse_lbfs = 9000; % lbf*s
impulse = impulse_lbfs * 4.44822; % N*s

% Calculate burn time
t_b = impulse / thrust; % s

% Calculate propellant mass
m_prop_0 = impulse / (Isp * g_0); % kg

% Calculate mass flow rate
m_dot = m_prop_0 / t_b; % kg/s

% Print title
fprintf('\n------ 1DOF SIMULATION ------\n\n')

% Print input parameters
fprintf('INPUT PARAMETERS\n')
fprintf('Dry Mass = %.2f lbm\n', m_dry_lbm)
fprintf('Thrust   = %.1f lbf\n', thrust_lbf)
fprintf('Isp      = %.1f s\n', Isp)
fprintf('Impulse  = %.1f lbf*s\n', impulse_lbfs)
fprintf('Cd       = %.2f\n', Cd)
fprintf('S        = %.4f m^2\n\n', S)

% Print calculated parameters
fprintf('CALCULATED PARAMETERS\n')
fprintf('Propellant Mass = %.2f lbm\n', m_prop_0 * 2.20462)
fprintf('Total Mass      = %.2f lbm\n', (m_prop_0+m_dry)*2.20462)
fprintf('Burn Time       = %.2f s\n', t_b)
fprintf('Mass Flow Rate  = %.2f lbm/s\n\n', m_dot * 2.20462)

% Initial conditions
alt = 0; % m
v = 0; % m/s
acc_old = 0; % m/s^2

% Time intervals
t = 0; % s
dt = 0.001; % s

% Maximum values
v_max = 0; % m/s
Mach_num_max = 0;
acc_max = 0; % m/s^2

fprintf('Simulation running...')
tic
while v >= 0
    
    % Increment time
    t = t + dt; % s
    
    % Calculate total mass
    if t < t_b
        m_prop = m_prop_0 - m_dot*t; % kg
    else
        m_prop = 0; % kg
    end
    m = m_dry + m_prop; % kg
    
    % Determine speed of sound and density at current altitude
    [~, a, ~, rho] = altcond1(alt);
    
    % Calculate drag force
    drag = 0.5 * rho * v^2 * Cd * S; % N
    
    % Calculate instantaneous acceleration
    if t < t_b
        acc_i = ((thrust - drag) / m) - g_0; % m/s^2
    else
        acc_i = (-drag / m) - g_0; % m/s^2
    end
    
    % Calculate second order acceleration
    acc = acc_i + (acc_i - acc_old) / 2; % m/s
    
    % Determine if max acceleration has been achieved
    if acc > acc_max
        acc_max = acc; % m/s^2
    end
        
    % Integrate acceleration
    v = v + acc*dt; % m/s
    
    % Assign current instantaneous acceleration to old variable
    acc_old = acc_i; % m/s
    
    % Determine if max speed has been achieved
    if v > v_max
        v_max = v; % m/s
        v_max_alt = alt; % m
    end
    
    % Calculate Mach number
    Mach_num = v / a;
    
    % Determine if maximum Mach number has been achieved
    if Mach_num > Mach_num_max
        Mach_num_max = Mach_num;
        Mach_num_max_alt = alt; % m
    end
    
    % Integrate velocity
    alt = alt + v*dt;
    
end
run_time = toc;
fprintf('done.\n')
fprintf('Simulation time = %.3f s\n\n', run_time)

% Print results
fprintf('RESULTS\n')
fprintf('Max Altitude     = %.2f ft\n', alt)
fprintf('Time to Apogee   = %.1f s\n', t);
fprintf('Max Speed        = %.1f m/s\n', v_max)
fprintf('Altitude at Max Speed = %.2f ft\n', v_max_alt/0.3048)
fprintf('Max Mach number  = %.2f \n', Mach_num_max)
fprintf('Altitude at Max Speed = %.2f ft\n', Mach_num_max_alt/0.3048)
fprintf('Max Acceleration = %.2f g''s \n\n', acc_max/g_0)