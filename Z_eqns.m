function [xdot] = Z_eqns(t,x,w)
global V X Y phiY m J cu Ma ts S T rho UY CG

Z = x(1);
phiZ = x(2);
UZ = w;

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
v_Z = [v_body(1); v_body(2); 0];
l_body = B2E\l;
l_Z = [l_body(1); l_body(2); 0];
d_body = B2E\d;
d_Z = [d_body(1); d_body(2); 0];
u_Z = [0; UZ; 0];

Z_dot = 1/J(3,3) * (cross(l_Z,cp) + cross(d_Z,cp) + cross(u_Z,cu));
Z_dot = Z_dot(3);

% phi angles in planes perpendicular to Y and Z body
phiZ_prev = phiZ;
phiZ = atan(v_Z(2)/v_Z(1)) + Z;
phiZ_dot = (phiZ - phiZ_prev)/ts;

xdot(1) = Z_dot;
xdot(2) = phiZ_dot;

end