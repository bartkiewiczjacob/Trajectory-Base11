function [xdot] = pitch_sim_eqns(t,x,w)

global T D L G V m J h cp cu ts

phi = x(1);
theta = x(2);
U = w;

B2E = [cos(pi/2+theta) 0 sin(pi/2+theta); 0 1 0; -sin(pi/2+theta) 0 cos(pi/2+theta)];
V2E = [cos(pi/2+phi) 0 sin(pi/2+phi); 0 1 0; -sin(pi/2+phi) 0 cos(pi/2+phi)];

att = B2E*[1; 0; 0];
v = V*V2E*[1; 0; 0]; 
v_unit = v/V;
v_perp = [-v(3); 0; v(1)];

t = T * att;
l = L * v_perp;
d = D * -v_unit;
g = [0; 0; G];
u = B2E*[0; 0; U];

f = t+l+d+g+u;
a = f/m;
h = a(3)*ts^2/2 + v(3)*ts + h;

v = v + a*ts;
V = norm(v);
phi_dot = (atan(v(3)/v(1)) - phi)/ts; % phi_dot = cos(phi)^2*[a(1)*v(3)-a(3)*v(1)]/v(3)^2;
theta_dot = 1/J(2,2) * (cross(l,cp) + cross(d,cp) + cross(u,cu));
theta_dot = theta_dot(2);

xdot(1) = theta_dot;
xdot(2) = phi_dot;

end
