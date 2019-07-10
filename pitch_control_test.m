thrust = 5000; % N
pitch0 = 3*pi/180; % rad
h0 = 0; % m
m_dot = -20*0.453592; % kg/s
diam = 18*0.0254; % m
r = diam/2; % m
len = 24*0.3048; % m
burn = 200000/T; % s
m_wet = 1300*0.453592; % kg
m_dry = m_wet + m_dot*burn; % kg
J_wet = m_wet * [1/12*len^2 + 1/4*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/2*r^2]; % kg m2
J_dry = m_dry * [ /12*len^2 + 1/4*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/2*r^2]; % kg m2
 
t_step = 0.5; % s
t0 = 0; % s
m = m_wet;
pitch = pitch0;
att = [sin(pitch) -cos(pitch)];
V = [0 0]; 
V_norm = V;
inc = acos(dot(V_norm, att));

while t <= burn
    T = thrust*att; % N, vector
    D = 