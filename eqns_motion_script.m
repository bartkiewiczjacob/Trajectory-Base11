clear; clc;

% retrieve aerodynamic data
Cd_data = load('cd_M.csv');
Cn_data = load('cl_M_a.csv');
Cp_data = load('cp_M_a.csv');
alpha_data = load('alpha_data.csv');
M_data = load('M_data.csv');

% design parameters
T = 5100*4.44822; % N
d = 18*0.0254; % m
r = d/2; % m
S = pi*d^2/4; % m2
len = 24*0.3048; % m

burn_time = 40; % s
m_wet = 1444.12*0.453592; % kg
mdot = 23.35*0.453592; % kg/s
m_dry = m_wet - mdot*burn_time; % kg
CG = len/2; % m

% initial position
theta0 = 0; % rad
psi0 = 0; % rad
h0 = 4595*0.3048; % m
u0 = 0; % m/s
v0 = 0;
w0 = 0;

% wind
lat = 32.9861;
lon = -106.9717;
day = 119;
sec = 12*3600;
h_trailer = 75*0.3048;

sim('eqns_motion_mdl.slx')

