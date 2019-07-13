clear; clc; close all;

global T D L G V m J h cp cu init ts

T = 10000; % N
m_dot = -20*0.453592; % kg/s
diam = 18*0.0254; % m
S = pi*diam^2/4; % m2
r = diam/2; % m
len = 24*0.3048; % m
burn = 200000/T; % s
m_wet = 1300*0.453592; % kg
m_dry = m_wet + m_dot*burn; % kg
J_wet = m_wet * [1/12*len^2 + 1/4*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/2*r^2]; % kg m2
J_dry = m_dry * [1/12*len^2 + 1/4*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/2*r^2]; % kg m2
CG_dot = -0.0557743*0.3048; % m 
CG_lim = len - 14.47*0.3048; % m
CU = 2; % m from tip
alpha_e = 0;
phi_e = 0;
u_e = 0;

global Ma_val CD_data_mach CL_data_mach CP_data_mach CD_fit_alpha CL_fit_alpha CP_fit_alpha
data_mach = load('aero_data_mach.csv');
data_alpha = load('aero_data_incidence.csv');
Ma_val = [0:0.01:5]';
alpha_val = [0:15]'; % rad
CD_data_mach = data_mach(:, 3);
CL_data_mach = data_mach(:, 8);
CP_data_mach = data_mach(:, 13)*0.0254;
CD_fit_alpha = [];
CL_fit_alpha = [];
CP_fit_alpha = [];
for iter = 1:5
    CD_fit_alpha = [CD_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,3), 2)'];
    CL_fit_alpha = [CL_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,8), 1)'];
    CP_fit_alpha = [CP_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,13)*0.0254, 4)'];
end
CD_fit_alpha(end, :) = [];
CL_fit_alpha(end, :) = [];
CP_fit_alpha(end, :) = [];

ts = 0.5; % s
h = 0; % mm = m_wet; % kg
m = m_wet;
J = m * [1/2*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/12*len^2 + 1/4*r^2];

theta = 0.1*pi/180; % rad
B2E = [cos(pi/2+theta) 0 sin(pi/2+theta); 0 1 0; -sin(pi/2+theta) 0 cos(pi/2+theta)];
phi = 0;
V = 0;
alpha = theta - phi;

CG = len - 16.6995*0.3048; % m from tip
tip = B2E*[-CG; 0; 0];
cp = B2E*[0; 0; 0]; 
CU = 1;
cu = B2E*[-(CG-CU); 0; 0];

U = 0;

figure(1); hold on;
time = 0;
plot(time, theta*180/pi, '*');
xlabel('Time [s]'); ylabel('Pitch [°]');
figure(2); hold on;
plot(time, V, '*')
xlabel('Time [s]'); ylabel('Ma [ ]');

%% first iteration
att = B2E*[1; 0; 0];

[temp, c, P, rho] = atmosisa(h);

Ma = round(V/c, 2);
[CD, CL, CP] = aero_coeff(Ma, alpha);
L = -sign(alpha)*1/2*rho*V^2*S*CL;
D = 1/2*rho*V^2*S*CD;
G = m*9.81;

v = [0; 0; 0];
v_unit = [0; 0; 0];
v_perp = [0; 0; 0];

t = T * att;
l = L * v_perp;
d = D * -v_unit;
g = [0; 0; G];
u = B2E*[0; 0; U];

f = t+l+d+g+u;
a = f/m;
h = a(3)*ts^2/2 + v(3)*ts + h;

theta_dot = 1/J(2,2) * (cross(l,cp) + cross(d,cp) + cross(u,cu));
theta_dot = theta_dot(2);

v = v + a*ts;
V = norm(v);

theta = theta + theta_dot*ts;
phi = atan(v(1)/v(3));
alpha = theta - phi;
init = [theta; phi]; % for sfun

m = m + m_dot*ts;
CG = CG + CG_dot*ts;
tip = B2E*[-CG; 0; 0];
cp = B2E*[CP-CG; 0; 0];
cu = B2E*[-(CG-CU); 0; 0];

[A,B,C,D]=linmod('pitch_sim_mdl',[0;0],0);
K = place(A,B,[-10 -11]);
U = -K*[theta; phi];

time = time + ts;
figure(1)
plot(time, theta*180/pi, '*');
figure(2)
plot(time, Ma, '*');

%% main loop
while time <= 20 
    [temp, c, P, rho] = atmosisa(h);

    Ma = round(V/c, 2);
    [CD, CL, CP] = aero_coeff(Ma, alpha);
    L = -sign(alpha)*1/2*rho*V^2*S*CL;
    D = 1/2*rho*V^2*S*CD;
    G = m*9.81;
    
    x(1) = theta;
    x(2) = phi;
    [xdot] = pitch_sim_eqns([],x,U);
    theta_dot = xdot(1);
    phi_dot = xdot(2);
    
    theta = theta + theta_dot*ts;
    phi = phi + phi_dot*ts;
    alpha = theta - phi;
    init = [theta; phi]; % for sfun
    
    m = m + m_dot*ts;
    CG = CG + CG_dot*ts;
    tip = B2E*[-CG; 0; 0];
    cp = B2E*[CP-CG; 0; 0];
    cu = B2E*[-(CG-CU); 0; 0];
    
    [A,B,C,D]=linmod('pitch_sim_mdl',[0;0],0);
    K = place(A,B,[-10 -11]);
    U = -K*[theta; phi];

    time = time + ts;
    if time == 2.5
        keyboard;
    end
    figure(1)
    plot(time, theta*180/pi, '*');
    figure(2)
    plot(time, Ma, '*');

end