% AUTHOR: Andy Meyer
% LAST MODIFIED: 11/07/2019

% FUNCTION PURPOSE
% Compute mass properties based on the time from t-zero

% INPUTS
% t - time in seconds from t-zero

% OUTPUTS
% mass   - mass in kilograms
% CG     - location of center of gravity in body frame in meters
%          [xb, yb, zb]
% I_vect - moments and products of inertia in kg*m^2
%          [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]

function [mass, CG, I_vect] = mass_func(t)
    
    % Initial propellant mass
    m_prop_0 = 21; % kg
    
    % Dry mass
    m_dry = 36; % kg

    % Burn time
    t_b = 9; % s
    
    % Mass flowrate
    m_dot = m_prop_0 / t_b; % kg/s

    % Calculate amount of propellant mass left
    if t < t_b
        m_prop = m_prop_0 - m_dot*t; % kg
    else
        m_prop = 0; % kg
    end
    
    % Determine total mass
    mass = m_dry + m_prop; % kg
    
    % Center of gravity
    CG = [0.6, 0, 0]; % m
    
    % Moment and product of inertia vector
    I_vect = [0.1372, 88, 88, 0, 0, 0]; % kg*m^2
    
end