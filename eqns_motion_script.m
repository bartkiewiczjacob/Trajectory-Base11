clear; clc;

% retrieve aerodynamic data
load('aerodata.mat');
air_data = load('air_data.csv');

% design parameters
d = 16*0.0254; % m
r = d/2; % m
S = pi*d^2/4; % m2
len = 24*0.3048; % m
throttle = 1;
T = 5100; %[lbs]
T_th = T*throttle;
t_th = 0;
margin = 60;
burn_time = (200000 - T*t_th)/T_th; % s
m_wet = 1437.33*0.453592 + margin*0.453592; % kg
mdot = 23.34*0.453592; % kg/s
mdot_th = mdot*throttle;
m_dry = m_wet - mdot*t_th - mdot_th*burn_time + margin*0.453592; % kg
CG = len/2; % m

% initial position
theta0 = 0; % rad
psi0 = 0; % rad
h0 = 4595*0.3048; % m
u0 = 0.1; % m/s
v0 = 0.0;
w0 = 0.0;

% wind
lat = 32.9861;
lon = -106.9717;
day = 119;
sec = 12*3600;
h_trailer = 75*0.3048;

%rail
lr = 75*0.3048;
x1 = len*3/4;
x2 = len*1/4;
mu = 0.42;

% engine %[lbs]
De = 8.4932;
Dt = 3.0968;
At = pi*Dt^2/4;
Pc = 450;
Pe = 10;
AR = De^2/Dt^2;

%heating data
cp = 910; %Material Specific Heat [J/kg*K]
material_density = 2700; %kg/m^3
x = 0.1; %Analysis Location from Nose Tip [m]
t = 0.07*0.0254; %Skin thickness[m]
kc = 0.35; %Insulation Thermal Conductivity [W/mK]
tc = 0.000; %Insulation Thickness [m]
k_he = 0.1177; %Helium Thermal Conductivity [W/mK]


 

