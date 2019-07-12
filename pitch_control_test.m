T = 5000; % N
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
CP_data_mach = data_mach(:, 13);
CD_fit_alpha = [];
CL_fit_alpha = [];
CP_fit_alpha = [];
for iter = 1:5
    CD_fit_alpha = [CD_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,3), 2)']
    CL_fit_alpha = [CL_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,8), 1)']
    CP_fit_alpha = [CP_fit_alpha, polyfit(alpha_val, data_alpha(iter:5:end,13), 4)']
end
CD_fit_alpha(end, :) = []
CL_fit_alpha(end, :) = []
CP_fit_alpha(end, :) = []

t_step = 0.5; % s
t = 0; % s
h = 0; % m
m = m_wet; % kg
v = [0 0]; % m/s
pitch = 3*pi/180; % rad
phi = 0;
alpha = pitch;
att = [sin(pitch) -cos(pitch)];
alpha_dot = 0; % rad/s
phi_dot = 0; % rad/s
CG = len - 16.6995*0.3048; % m from tip
CP = CG; 
l = CG - CU; % m
ksi = CP - CG; % m
U = 0;

% first iteration
V = norm(v);
v_perp = [v(2), -v(1)];
v_unit = v/V;

[T, c, P, rho] = atmosisa(h);
Ma = V/c;
[CD, CL, CP] = aero_coeff(Ma, alpha);

t = T*att;
L = 1/2*rho*V^2*S*CL;
l = L * v_perp;
D = 1/2*rho*V^2*S*CD;
d = D * -v_unit;
g = [0, m*g];

f = t+l+d+g;
a = f/m;
h = a*t_step^2/2 + V*t_step + h;
v = v + a*t_step;

alpha_dot = 1/J(2,2) * (-L*cos(alpha)*ksi + D*sin(alpha)*ksi + U*l);
alpha = alpha + alpha_dot*t_step;

m = m + m_dot*t_step;
J = m * [1/2*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/12*len^2 + 1/4*r^2];
CG = CG + CG_dot*t_step;
l = CG - CU; % m
ksi = CP - CG; % m

while t <= burn
    [T, c, P, rho] = atmosisa(h); 
    
    if m + m_dot*t_step >= m_dry
        m = m + m_dot*t_step;
    else
        m = m_dry;
    end
    
    if CG + CG_dot*t_step >= CG_lim
        CG = CG + CG_dot*t_step;
    else
        CG = CG_lim;
    end
        
    J = m * [1/2*r^2, 0, 0; 0, 1/12*len^2 + 1/4*r^2, 0; 0, 0, 1/12*len^2 + 1/4*r^2];
    
    Ma = norm(V)/c;
    
    
    T = thrust*att; % N, vector
    