

% Script initialization
clear, clc

% Intertial Frame
% [X, Y, Z] = [East, North, Up]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravitational acceleration
g = 9.8; % m/s^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temporary spot for aerodynamic properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_ax = 0;
C_ny = 0;
C_nz = 0;
Sd = 0;
Sl = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define initial orientation vectors of rocket in inertial frame
xb_hat_0 = [1, 0, 0];
yb_hat_0 = [0, 1, 0];
zb_hat_0 = [0, 0, 1];

% Initial quaternion
q_0 = [1, 0, 0, 0];
q = q_0;

% Rotate rocket to up-east-north
theta = deg2rad(-90); % rad
rot_vect = [0, 1, 0]; % intertial y-axis
q_v = (rot_vect / norm(rot_vect)) .* sin(theta/2);
q_s = cos(theta/2);
q = Quat_Mult([q_s, q_v], q);
theta = deg2rad(-90); % rad
rot_vect = [0, 0, 1]; % inertial z-axis
q_v = (rot_vect / norm(rot_vect)) .* sin(theta/2);
q_s = cos(theta/2);
q = Quat_Mult([q_s, q_v], q);

% Azimuthal angle
azimuth_deg = 0; % deg
azimuth = azimuth_deg * pi/180; % radians
theta = -azimuth; % rad
rot_vect = [0, 0, 1];
q_v = (rot_vect / norm(rot_vect)) .* sin(theta/2);
q_s = cos(theta/2);
q = Quat_Mult([q_s, q_v], q);

% Declination angle
declination_deg = 90; % deg
theta = -deg2rad(90 - declination_deg); % rad
rot_vect = Quat_Rot(yb_hat_0, q);
q_v = (rot_vect / norm(rot_vect)) .* sin(theta/2);
q_s = cos(theta/2);
q = Quat_Mult([q_s, q_v], q);

% Roll angle
roll_angle_deg = 0; % deg
theta = deg2rad(roll_angle_deg); % rad
rot_vect = Quat_Rot(xb_hat_0, q);
q_v = (rot_vect / norm(rot_vect)) .* sin(theta/2);
q_s = cos(theta/2);
q = Quat_Mult([q_s, q_v], q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define time interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time Interval
t_0 = 0; % seconds
t_f = 110; % seconds
dt = 0.0001; % seconds
time_list = t_0:dt:t_f; % seconds
numTimeSteps = length(time_list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial latitude, longitude, and elevation
latitude = 40.5; % degrees
longitude = -87; % degrees
elevation = 753*0.3048; % m

% Initial position
p_x0 = 0; % m
p_y0 = 0; % m
p_z0 = 0; % m

% Initial velocity
v_x0 = 0; % m/s
v_y0 = 0; % m/s
v_z0 = 0; % m/s

% Initial angular velocity in body frame
w_b_x0 = 0; % rad/s
w_b_y0 = 0; % rad/s
w_b_z0 = 0; % rad/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Position variable
p_x = p_x0; % m
p_y = p_y0; % m
p_z = p_z0; % m
p = [p_x, p_y, p_z]; % [m, m, m]

% Velocity variable
v_x = v_x0; % m/s
v_y = v_y0; % m/s
v_z = v_z0; % m/s
v = [v_x, v_y, v_z]; % [m/s, m/s, m/s]

% Angular velocity variable
w_x = w_b_x0; % rad/s
w_y = w_b_y0; % rad/s
w_z = w_b_z0; % rad/s
w_b = [w_x, w_y, w_z]; % [rad/s, rad/s, rad/s]

% Angular acceleration
w_dot_x_old = 0; % rad/s^2
w_dot_y_old = 0; % rad/s^2
w_dot_z_old = 0; % rad/s^2

% Acceleration
a_x_old = 0; % m/s^2
a_y_old = 0; % m/s^2
a_z_old = 0; % m/s^2

% Arrays
p_list = zeros(round(numTimeSteps/50), 3);

data_index = 0;
for index = 1:numTimeSteps
    t = time_list(index); % s
    
    % Check to ensure altitude > 0
    if p_z < 0
        fprintf('Altitude < 0 m\n')
        break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mass Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Call function to determine mass properties
    [mass, CG_b, I_vect] = mass_func(t); % kg
    Ixx = I_vect(1); % kg*m^2
    Iyy = I_vect(2); % kg*m^2
    Izz = I_vect(3); % kg*m^2
        % mass - kilograms at current time
        % CG_b - center of gravity location [xb, yb, zb] in meters
        % I_vect - moments and products of inertia about body frame
        %          axes, [Ixx, Iyy, Izz, Ixy, Ixz, Iyz], kg*m^2
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Thrust Forces and Moments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Call function to determine current thrust
    thrust = thrust_func(t); % N
    
    % Assign point through which thrust acts
    xb_thrust = 0.1; % m
    yb_thrust = 0.0001; % m
    zb_thrust = 0; % m
    p_thrust = [xb_thrust, yb_thrust, zb_thrust];
    
    % Assign direction in which thrust is acting
    xb_d_thrust = 1;
    yb_d_thrust = 0;
    zb_d_thrust = 0;
    d_thrust = [xb_d_thrust,yb_d_thrust,zb_d_thrust] ...
                 / norm([xb_d_thrust,yb_d_thrust,zb_d_thrust]);
             
	% Compute thrust components in rocket axes
    Thrust_b = thrust * d_thrust; % [N, N, N]
    
    % Compute position vector of thrust relative to CG
    r_thrust = p_thrust - CG_b; % [m, m, m]
    
    % Compute moments in rocket axes due to thrust
    M_thrust = cross(r_thrust, Thrust_b); % [N*m, N*m, N*m]
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aerodynamic Forces and Moments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Call function to determine atmospheric properties
    altitude = p_z; % m
    alt_geopotential = Geometric_to_Geopotential(altitude); % m
    [~, a_atm, ~, rho_atm] = atmosisa(alt_geopotential);
        % a_atm - speed of sound, m/s
        % rho_atm - atmospheric density, kg/m^3
        
	% Define wind velocity in inertial frame
    wind_i = [0, 0, 0]; % [m/s, m/s, m/s]
    
    % Wind relative rocket velocity
    v_wr = v - wind_i; % [m/s, m/s, m/s]
    v_inf = norm(v_wr); % m/s
    
    % Dynamic pressure
    Q_inf = 0.5*rho_atm*v_inf^2;

    % Wind velocity in rocket frame
    V_w_b = Quat_Rot(-v_wr, Quat_Inv(q));

    % Rocket axial orientation vector
    xb_hat = Quat_Rot(xb_hat_0, q);
    
    % Absolute angle of attack
    if norm(v_wr)
        AoA_deg = acosd(dot(xb_hat,v_wr)/(norm(xb_hat)*norm(v_wr)));
    else
        AoA_deg = 0;
    end

    % Pitch angle (angle of attack) and yaw angle (side-slip angle)
    V_w_bx = V_w_b(1); % m/s
    V_w_by = V_w_b(2); % m/s
    V_w_bz = V_w_b(3); % m/s
    if V_w_bx == 0
        alpha_deg = 0;
        beta_deg = 0;
    else
        alpha_deg =  atand(V_w_bz/V_w_bx); % deg
        beta_deg  = -atand(V_w_by/V_w_bx); % deg
    end
    
    % Drag
    drag = Q_inf * 0.4 * pi*(3*0.0254)^2; % N
    if drag ~= 0
        Drag_X = drag * -(v_wr(1)/v_inf); % N
        Drag_Y = drag * -(v_wr(2)/v_inf); % N
        Drag_Z = drag * -(v_wr(3)/v_inf); % N
    else
        Drag_X = 0;
        Drag_Y = 0;
        Drag_Z = 0;
    end

    % Call function to deterine aerodynamic forces
    %[A_x, N_y, N_z, CP_x, CP_y, CP_z] = ...
    %        aero_func(AoA_deg, alpha_deg, beta_deg, Q_inf, a_atm);
    A_x =  Q_inf*v_inf^2 * C_ax * Sd; % N
    N_y =  Q_inf * C_ny * Sl; % N
    N_z = -Q_inf * C_nz * Sl; % N
    CP_x = CG_b(1)-0.15; % m
    CP_y = CG_b(1)-0.15; % m
    CP_z = CG_b(1)-0.15; % m
    
    % Compute position vectors of CP's relative to CG
    r_Ax = CP_x - CG_b; % [m, m, m]
    r_Ny = CP_y - CG_b; % [m, m, m]
    r_Nz = CP_z - CG_b; % [m, m, m]
    
    % Compute aerodynamic moments
    M_aero_x = cross(r_Ax, [A_x, 0, 0]); % N*m
    M_aero_y = cross(r_Ny, [0, N_y, 0]); % N*m
    M_aero_z = cross(r_Nz, [0, 0, N_z]); % N*m
    M_aero = M_aero_x + M_aero_y + M_aero_z; % [N*m,N*m,N*m]
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Angular Velocity and Quaternion Integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find sum of moments about CG
    M_total = M_thrust + M_aero; % [N*m, N*m, N*m]
    M_x = M_total(1); % N*m
    M_y = M_total(2); % N*m
    M_z = M_total(3); % N*m
    
    % Compute instantaneous angular accelerations
    w_dot_x = M_x / Ixx; % rad/s^2
    w_dot_y = M_y / Iyy; % rad/s^2
    w_dot_z = M_z / Izz; % rad/s^2
    
    % Compute 2nd order angular accelerations over time step
    w_dot_x_avg = w_dot_x + (w_dot_x-w_dot_x_old)/2; % rad/s^2
    w_dot_y_avg = w_dot_y + (w_dot_y-w_dot_y_old)/2; % rad/s^2
    w_dot_z_avg = w_dot_z + (w_dot_z-w_dot_z_old)/2; % rad/s^2
    
    % Assign current angular acceleration to old variable
    w_dot_x_old = w_dot_x; % rad/s^2
    w_dot_y_old = w_dot_y; % rad/s^2
    w_dot_z_old = w_dot_z; % rad/s^2    
    
    % Integrate angular accelerations over time step
    w_x = w_x + w_dot_x_avg*dt; % rad/s
    w_y = w_y + w_dot_y_avg*dt; % rad/s
    w_z = w_z + w_dot_z_avg*dt; % rad/s
    w_b = [w_dot_x_avg, w_dot_y_avg, w_dot_z_avg]; % [rad/s]
    
    % Integrate quaternion over time step
    q_dot = 0.5 * Quat_Mult(q,[0,w_b]);
    q_old = q;
    q = q + q_dot*dt;
    q = q ./ norm(q);
    
    % Update rocket frame location (p) based on new quaternion
    p = Quat_Rot(CG_b,q_old) - Quat_Rot(CG_b,q) + p; % [m, m, m]
    p_x = p(1); % m
    p_y = p(2); % m
    p_z = p(3); % m
    
    % Check to ensure altitude > 0
    if p_z < 0
        fprintf('Altitude < 0 m\n')
        break
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocity and Position Integration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Rotate thrust force into inertial frame
    Thrust_i = Quat_Rot(Thrust_b, q); % [N, N, N]
    
    % Rotate aerodynamic forces into inertial frame
    % Aero_i = Quat_Rot([A_x, N_y, N_z], q); % [N, N, N]
    
    % Compute acceleration in each inertial direction
    a_x = (Thrust_i(1) + Drag_X)/mass; % m/s^2
    a_y = (Thrust_i(2) + Drag_Y)/mass; % m/s^2
    a_z = (Thrust_i(3) + Drag_Z)/mass - g; % m/s^2
    
    % Compute 2nd order acceleration over time step
    a_x_avg = a_x + (a_x-a_x_old)/2; % m/s^2
    a_y_avg = a_y + (a_y-a_y_old)/2; % m/s^2
    a_z_avg = a_z + (a_z-a_z_old)/2; % m/s^2
    
    % Assign current accleration to old variable
    a_x_old = a_x; % m/s^2
    a_y_old = a_y; % m/s^2
    a_z_old = a_z; % m/s^2
    
    % Integrate acceleration to find velocity
    v_x = v_x + a_x_avg*dt; % m/s
    v_y = v_y + a_y_avg*dt; % m/s
    v_z = v_z + a_z_avg*dt; % m/s
    v = [v_x, v_y, v_z]; % [m/s, m/s, m/s]
    
    % Integrate velocity to find position
    p_x = p_x + v_x*dt; % m
    p_y = p_y + v_y*dt; % m
    p_z = p_z + v_z*dt; % m
    p = [p_x, p_y, p_z]; % [m, m, m]
    if mod(index,100) == 0
        data_index = data_index + 1;
        p_list(data_index, 1:3) = [p_x, p_y, p_z]; % [m, m, m]
    end
    
end

% Remove empty rows from data matrices
p_list = p_list(1:data_index,:);

% Convert p_list to latitude, longitude, and distance above ground
LLA = zeros(data_index,3);
LLA(:,1) = (p_list(:,2) ./ 111194.9) + latitude; % degrees
LLA(:,2) = (p_list(:,1) ./ 111194.9) + longitude; % degrees
LLA(:,3) = p_list(:,3) + elevation; % m


% Compute rocket orientation vectors at final time step
xb_hat = Quat_Rot(xb_hat_0, q);
yb_hat = Quat_Rot(yb_hat_0, q);
zb_hat = Quat_Rot(zb_hat_0, q);

% Plot orientation of rocket at final time step
figure(1)
clf
hold on
plot3([0,1.5*xb_hat(1)], [0,1.5*xb_hat(2)], [0,1.5*xb_hat(3)], 'r')
plot3([0,yb_hat(1)], [0,yb_hat(2)], [0,yb_hat(3)], 'g')
plot3([0,zb_hat(1)], [0,zb_hat(2)], [0,zb_hat(3)], 'b')
xlabel('X (East)')
ylabel('Y (North)')
zlabel('Z (Up)')
xlim([-2, 2])
ylim([-2, 2])
zlim([-2, 2])
grid on
daspect([1, 1, 1])
view([30, 34])
hold off

% Create Google Earth trajectory file
GreenTrajectoryPlotter(LLA(:,1),LLA(:,2),LLA(:,3));



