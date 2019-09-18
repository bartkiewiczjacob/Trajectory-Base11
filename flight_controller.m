clear; close all;

% make the following variables accessible to all scripts and functions
% without needing to define them
global V X Y Z phiY phiZ m J cu Ma ts S T rho UY UZ CG initY initZ time

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

% body to Earth initial error angles
X = 0; % roll angle
Y = 0.1*pi/180; % pitch angle
Z = Y; % yaw angle

% B2E matrix
% Earth frame: x left, y in, z down
% Body frame (initially): x up, y in, z left
BX = [1, 0, 0; 0, cos(X), -sin(X); 0, sin(X), cos(X)];
BY = [cos(pi/2+Y), 0, sin(pi/2+Y); 0, 1, 0; -sin(pi/2+Y), 0, cos(pi/2+Y)];
BZ= [cos(Z), -sin(Z), 0; sin(Z), cos(Z), 0; 0, 0, 1];
B2E = BX*BY*BZ; % order of rotations: X, Y, Z

att = B2E*[1; 0; 0]; % attitude vector
v = [0; 0; 0];
V = norm(v);
v_unit = v;
v_perp = v;

% distances in body frame
CG = len - 16.6995*0.3048; % tip to CG, m
tip = [-CG; 0; 0]; % tip to CG, vector
CP = 0; % m
cp = [CP-CG; 0; 0]; % CP to CG, vector
CU = 1; % tip to CU, m
cu = [-(CG-CU); 0; 0]; % CU to CG, vector

UY = 0; % input magnitude, N
UZ = 0;
L = 0;
D = 0;
G = m*9.81;

time = 0;

% initialize plots
figure(1); hold on;
plot(time, Y*180/pi, '*g', time, Z*180/pi, '*r');
xlabel('Time [s]'); ylabel('Angle [°]');
% figure(2); hold on;
% plot(time, Z*180/pi, '*b')
% xlabel('Time [s]'); ylabel('Z [°]');
% figure(3); hold on;
% xlabel('Time [s]'); ylabel('Ma [ ]');
% figure(4); hold on;
% xlabel('Time [s]'); ylabel('\Phi [°]');
figure(5); hold on;
xlabel('Time [s]'); ylabel('Control moment [Nm]');

%% main loop
while time <= 28 
    % forces 
    t = T * att;
    l = L * v_perp;
    d = D * -v_unit;
    g = [0; 0; G];
    u = B2E*[0; UZ; UY];
    
    % kinematics, assume no dynamics during ts
    f = t+l+d+g+u;
    a = f/m;
    h = -a(3)*ts^2/2 - v(3)*ts + h;
    v = v+a*ts;
    V = norm(v);
    v_unit = v/V;
    h = -a(3)*ts^2/2 - v(3)*ts + h;
    
    % projections in planes perpendicular to pitch (Y) and yaw (Z) axes
    v_body = B2E\v;
    v_Y = [v_body(1); 0; v_body(3)];
    v_Z = [v_body(1); v_body(2); 0];
    l_body = B2E\l;
    l_Y = [l_body(1); 0; l_body(3)];
    l_Z = [l_body(1); l_body(2); 0];
    d_body = B2E\d;
    d_Y = [d_body(1); 0; d_body(3)];
    d_Z = [d_body(1); d_body(2); 0];
    u_Y = [0; 0; UY];
    u_Z = [0; UZ; 0];
    
    % angular rates
    X_dot = 0;
    X = X + X_dot*ts;
    Y_dot = 1/J(2,2) * (cross(l_Y,cp) + cross(d_Y,cp) + cross(u_Y,cu));
    Y_dot = Y_dot(2);
    Y = (Y + Y_dot*ts);
    Z_dot = 1/J(3,3) * (cross(l_Z,cp) + cross(d_Z,cp) + cross(u_Z,cu));
    Z_dot = Z_dot(3);
    Z = (Z + Z_dot*ts);
    
    % angles between velocity and vertical in planes perpendicular to Y and Z
    phiY_prev = phiY;
    phiZ_prev = phiZ;
    v_Ye = B2E*v_Y;
    v_Ze = B2E*v_Z;
    phiY = atan(v_Ye(1)/v_Ye(3));
    phiZ = atan(-v_Ze(2)/v_Ze(3));
    phiY_dot = (phiY - phiY_prev)/ts;
    phiZ_dot = (phiZ - phiZ_prev)/ts;

    % initial states for s-functions
    initY = [Y; phiY];
    initZ = [Z; phiZ];
    
    % B2E matrices
    BX = [1, 0, 0; 0, cos(X), -sin(X); 0, sin(X), cos(X)];
    BY = [cos(pi/2+Y), 0, sin(pi/2+Y); 0, 1, 0; -sin(pi/2+Y), 0, cos(pi/2+Y)];
    BZ= [cos(Z), -sin(Z), 0; sin(Z), cos(Z), 0; 0, 0, 1];
    B2E = BX*BY*BZ; % order of rotations: X, Y, Z
    
    att = B2E*[1; 0; 0]; % attitude vector
    alpha = acos(dot(att,v_unit));
    v_perp = cross(v_unit, cross(att,v_unit));
    v_perp = v_perp/norm(v_perp);
    
    % physical properties
    m = m + m_dot*ts;
    J = m * [1/2*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/12*len^2 + 1/4*r^2]; % kg m2
    
    [temp, c, P, rho] = atmosisa(h);
    Ma = round(V/c, 2);
    
    % aerodynamics
    [CD, CL, CP] = aero_coeff(Ma, abs(alpha));
    L = 1/2*rho*V^2*S*CL;
    D = 1/2*rho*V^2*S*CD;
    G = m*9.81;
    
    % CP and CU to CG vectors
    CG = CG + CG_dot*ts;
    tip = [-CG; 0; 0];
    cp = [CP-CG; 0; 0];
    cu = [-(CG-CU); 0; 0];

    [Ay,By,Cy,Dy]=linmod('Y_sfun_mdl', [theta_e;phi_e], u_e); % linearization of s-function describing [Y_dot, phiY_dot]=Y_eqns(Y, phiY)
    Ky = place(Ay, By, poles); % gains for desired poles of A-BK
    UY = -Ky*[Y; phiY]; % control input
    
    [Az,Bz,Cz,Dz]=linmod('Z_sfun_mdl', [theta_e;phi_e], u_e); % linearization
    Kz = place(Az, Bz, poles); % gains for desired poles of A-BK
    UZ = -Kz*[Z; phiZ]; % control input
    
    time = time + ts;
    
%     plots
    figure(1)
    plot(time, Y*180/pi, '*g', time, Z*180/pi, '*r');
%     figure(2)
%     plot(time, Z*180/pi, '*b');
%     figure(3); 
%     plot(time, Ma, '*b')
%     figure(4);
%     plot(time, phiY*180/pi, '*g', time, phiZ*180/pi, '*r')
    figure(5);
    plot(time, UY, '*g', time, UZ, '*r');
    
end

figure(1)
legend('Pitch', 'Yaw')
figure(5)
ylabel('Control force [N]')
legend('Pitch', 'Yaw')