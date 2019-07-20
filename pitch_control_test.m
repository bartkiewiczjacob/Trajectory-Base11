clear; clc; close all;

% make the following variables accessible to all scripts and functions
% without needing to define them
global T V m J cp cu init ts time Ma rho S CG

% rocket parameters and properties
T = 7000; % N
m_dot = -20*0.453592; % kg/s
diam = 18*0.0254; % m
S = pi*diam^2/4; % m2
r = diam/2; % m
len = 24*0.3048; % m
burn = 200000/T; % s
m_wet = 1300*0.453592; % kg
m_dry = m_wet + m_dot*burn; % kg
CG_dot = -0.0557743*0.3048; % m 
CG_lim = len - 14.47*0.3048; % m
CU = 2; % m from tip

% control parameters
theta_e = 0; % equilibrium points
phi_e = 0;
u_e = 0;
poles = [-2 -1]; % poles of A-BK

% retrieve data from aerodynamics test and make it easily accessible by
% other functions and scripts
global Ma_val CD_data_mach CL_data_mach CP_data_mach CD_fit_alpha CL_fit_alpha CP_fit_alpha
data_mach = load('aero_data_mach.csv');
data_alpha = load('aero_data_incidence.csv');
Ma_val = [0:0.01:5];
alpha_val = [0:15]'; % rad
% data varying with Mach number
CD_data_mach = data_mach(:, 3);
CL_data_mach = data_mach(:, 8);
CP_data_mach = data_mach(:, 13)*0.0254; 
% data varying with angle of attack is simpler so we get the polyfit for a
% few Mach numbers
CD_fit_alpha = [];
CL_fit_alpha = [];
CP_fit_alpha = [];
for iter = 1:5
    CD_fit_alpha = [CD_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,3), 2)'];
    CL_fit_alpha = [CL_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,8), 1)'];
    CP_fit_alpha = [CP_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,13)*0.0254, 4)'];
end
% remove the last coefficient of the polyfit because that one is set by the
% number at the specific Mach number and angle of attack = 0 (from
% C._data_mach)
CD_fit_alpha(end, :) = [];
CL_fit_alpha(end, :) = [];
CP_fit_alpha(end, :) = [];

%% initialize
ts = 0.5; % s
h = 0; % m
m = m_wet; % kg
J = m * [1/2*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/12*len^2 + 1/4*r^2]; % kg m2

theta = 0.1*pi/180; % pitch angle (attitude to vertical), rad
B2E = [cos(pi/2+theta) 0 sin(pi/2+theta); 0 1 0; -sin(pi/2+theta) 0 cos(pi/2+theta)]; % transformation matrix, body -> Earth frame
phi = 0; % velocity to vertical, rad
V = 0; % m/s
alpha = theta - phi; % angle of attack, rad
init = [theta; phi]; % initial state for s-function (not important)

CG = len - 16.6995*0.3048; % tip to CG, m
tip = B2E*[-CG; 0; 0]; % tip to CG, vector
CP = 0; % m
cp = B2E*[CP; 0; 0]; % CP to CG, vector
CU = 1; % tip to CU, m
cu = B2E*[-(CG-CU); 0; 0]; % CU to CG, vector

U = 0; % input magnitude, N
L = 0;
D = 0;
G = m*9.81;

time = 0;
U_rec = [U];

% initialize plots
figure(1); hold on;
plot(time, theta*180/pi, '*k');
xlabel('Time [s]'); ylabel('Pitch [°]');
figure(2); hold on;
plot(time, V, '*y')
xlabel('Time [s]'); ylabel('Ma [ ]');
figure(3); hold on;
plot(time, h, '*b')
xlabel('Time [s]'); ylabel('Altitude [m]');
figure(4); hold on;
xlabel('Time [s]'); ylabel('\Phi [°]');
figure(5); hold on;
xlabel('Time [s]'); ylabel('Control moment [Nm]');

%% main loop
while time <= burn 
    % transformation matrices
    B2E = [cos(pi/2+theta) 0 sin(pi/2+theta); 0 1 0; -sin(pi/2+theta) 0 cos(pi/2+theta)]; % body -> Earth
    V2E = [cos(pi/2+phi) 0 sin(pi/2+phi); 0 1 0; -sin(pi/2+phi) 0 cos(pi/2+phi)]; % velocity vector -> earth

    att = B2E*[1; 0; 0]; % attitude vector 
    v = V*V2E*[1; 0; 0]; % velocity vector
    % compute unit vectors
    if V~=0
        v_unit = v/V;
        v_perp = [-v(3); 0; v(1)]/V;
    else 
        v_unit = v;
        v_perp = v;
    end
    
    % forces 
    t = T * att;
    l = L * v_perp;
    d = D * -v_unit;
    g = [0; 0; G];
    u = B2E*[0; 0; U];
    
    % kinematics, assume no dynamics during ts
    f = t+l+d+g+u;
    a = f/m;
    h = -a(3)*ts^2/2 - v(3)*ts + h;
    v = v + a*ts;
    V = norm(v);
    
    % angles
    theta_dot = 1/J(2,2) * (cross(l,cp) + cross(d,cp) + cross(u,cu));
    theta_dot = theta_dot(2);
    theta = (theta + theta_dot*ts);
    phi = atan(v(1)/v(3));
    alpha = theta - phi;
    init = [theta; phi]; % initial state for s-function (not important)
    
    % physical properties
    m = m + m_dot*ts;
    J = m * [1/2*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/12*len^2 + 1/4*r^2]; % kg m2
    
    [temp, c, P, rho] = atmosisa(h);
    Ma = round(V/c, 2);
    
    % aerodynamics
    [CD, CL, CP] = aero_coeff(Ma, alpha);
    L = -sign(alpha)*1/2*rho*V^2*S*CL;
    D = 1/2*rho*V^2*S*CD;
    G = m*9.81;
    
    % CP and CU to CG vectors
    CG = CG + CG_dot*ts;
    tip = B2E*[-CG; 0; 0];
    cp = B2E*[CP-CG; 0; 0];
    cu = B2E*[-(CG-CU); 0; 0];
    
    [A,B,Cm,Dm]=linmod('pitch_sim_mdl', [theta_e;phi_e], u_e); % linearization
    K = place(A, B, poles); % gains for desired poles of A-BK
    U = -K*[theta; phi]; % control input
    U_rec = [U_rec; U];
    
    time = time + ts;
    
    % plots
    figure(1)
    plot(time, theta*180/pi, '*k');
    figure(2)
    plot(time, Ma, '*y');
    figure(3); 
    plot(time, h, '*b')
    figure(4);
    plot(time, phi*180/pi, '*r')
    figure(5);
    c_moment = cross(u,cu);
    plot(time, c_moment(2), '*g');

end