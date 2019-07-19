function [xdot] = pitch_sim_eqns(t,x,w)
% This function is called s-function pitch_sim_mdl (Simulink block
% written in code).
% This function is mainly the same as the main script, but it is in the
% form xdot = f(t,x,w), with:
% t = time
% x = state vector
% w = control input
% This is necessary in order to get the right state space matrices when we
% linearize

global T V m J cu ts Ma rho S CG% time varying parameters

%inputs
theta = x(1);
phi = x(2);
U = w;
alpha = theta - phi;

B2E = [cos(pi/2+theta) 0 sin(pi/2+theta); 0 1 0; -sin(pi/2+theta) 0 cos(pi/2+theta)];
V2E = [cos(pi/2+phi) 0 sin(pi/2+phi); 0 1 0; -sin(pi/2+phi) 0 cos(pi/2+phi)];

att = B2E*[1; 0; 0];
v = V*V2E*[1; 0; 0];
v_unit = v/V;
v_perp = [-v(3); 0; v(1)];

% aerodynamics
[CD, CL, CP] = aero_coeff(Ma, alpha);
L = -sign(alpha)*1/2*rho*V^2*S*CL;
D = 1/2*rho*V^2*S*CD;
G = m*9.81;
cp = B2E*[CP-CG; 0; 0];

t = T * att;
l = L * v_perp;
d = D * -v_unit;
g = [0; 0; G];
u = B2E*[0; 0; U];

f = t+l+d+g+u;
a = f/m;
v = v + a*ts;

phi_dot = (atan(v(1)/v(3)) - phi)/ts; % phi_dot = cos(phi)^2*[a(1)*v(3)-a(3)*v(1)]/v(3)^2;
theta_dot = 1/J(2,2) * (cross(l,cp) + cross(d,cp) + cross(u,cu));
theta_dot = theta_dot(2);

% outputs
xdot(1) = theta_dot;
xdot(2) = phi_dot;

end
