function [xdot] = Y_eqns(t,x,w)
global V X Z phiZ m J cu Ma ts S T rho UZ CG

Y = x(1);
phiY = x(2);
UY = w;

% B2E matrices
BX = [1, 0, 0; 0, cos(X), -sin(X); 0, sin(X), cos(X)];
BY = [cos(pi/2+Y), 0, sin(pi/2+Y); 0, 1, 0; -sin(pi/2+Y), 0, cos(pi/2+Y)];
BZ= [cos(Z), -sin(Z), 0; sin(Z), cos(Z), 0; 0, 0, 1];
B2E = BX*BY*BZ; % order of rotations: X, Y, Z
RX = [1, 0, 0; 0, cos(X), -sin(X); 0, sin(X), cos(X)];
RY = [cos(pi/2+phiY), 0, sin(pi/2+phiY); 0, 1, 0; -sin(pi/2+phiY), 0, cos(pi/2+phiY)];
RZ= [cos(phiZ), -sin(phiZ), 0; sin(phiZ), cos(phiZ), 0; 0, 0, 1];
V2E = RX*RY*RZ; % order of rotations: X, Y, Z

v = V*V2E*[1;0;0];
v_unit = v/V;

att = B2E*[1; 0; 0]; % attitude vector
alpha = acos( dot(att,v_unit) );
v_perp = cross(v_unit, cross(att,v_unit));
v_perp = v_perp/norm(v_perp);
    
% aerodynamics
[CD, CL, CP] = aero_coeff(Ma, abs(alpha));
L = 1/2*rho*V^2*S*CL;
D = 1/2*rho*V^2*S*CD;
G = m*9.81;
cp = [CP-CG; 0; 0]; % CP to CG, vector

% forces 
t = T * att;
l = L * v_perp;
d = D * -v_unit;
g = [0; 0; G];
u = B2E*[0; UZ; UY];

% kinematics, assume no dynamics during ts
f = t+l+d+g+u;
a = f/m;
v = v+a*ts;

% projections
v_body = B2E\v;
v_Y = [v_body(1); 0; v_body(3)];
l_body = B2E\l;
l_Y = [l_body(1); 0; l_body(3)];
d_body = B2E\d;
d_Y = [d_body(1); 0; d_body(3)];
u_Y = [0; 0; UY];

Y_dot = 1/J(2,2) * (cross(l_Y,cp) + cross(d_Y,cp) + cross(u_Y,cu));
Y_dot = Y_dot(2);

% phi angles in planes perpendicular to Y and Z body
phiY_prev = phiY;
phiY = atan(-v_Y(3)/v_Y(1)) + Y;
phiY_dot = (phiY - phiY_prev)/ts;

xdot(1) = Y_dot;
xdot(2) = phiY_dot;

end