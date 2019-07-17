function [xdot] = Z_eqns(t,x,w)

global V X Y phiY_body m J CG cu Ma UY ts rho S T

Z = x(1);
phiZ_body = x(2);
UZ = w;

% B2E matrices
BX = [1, 0, 0; 0, cos(X), -sin(X); 0, sin(X), cos(X)];
BY = [cos(pi/2+Y), 0, sin(pi/2+Y); 0, 1, 0; -sin(pi/2+Y), 0, cos(pi/2+Y)];
BZ= [cos(Z), -sin(Z), 0; sin(Z), cos(Z), 0; 0, 0, 1];
B2E = BX*BY*BZ; % order of rotations: X, Y, Z

Vx = sqrt( V/(1+tan(phiY_body)^2+tan(phiZ_body)^2) );
Vy = Vx*tan(Z);
Vz = Vx*tan(Y);

v_body = [Vx; Vy; Vz];
v = B2E*v_body;
v_unit = v/V;
att = B2E*[1; 0; 0];

theta = acos( dot(att, [-1; 0; 0]) );
phi = acos( dot(v_unit, [-1; 0; 0]) );
alpha = theta - phi;
v_perp = cross(att, cross(att,v_unit));
v_perp = v_perp/norm(v_perp);

% aerodynamics
[CD, CL, CP] = aero_coeff(Ma, abs(alpha));
L = 1/2*rho*V^2*S*CL;
D = 1/2*rho*V^2*S*CD;
G = m*9.81;
cp = [CP-CG; 0; 0];

% forces 
t = T * att;
l = L * v_perp;
d = D * -v_unit;
g = [0; 0; G];
u = B2E*[0; UZ; UY];

% kinematics, assume no dynamics during ts
f = t+l+d+g+u;
a = f/m;
v = v + a*ts;
v_unit = v/norm(v);

% projections
l_body = B2E\l;
l_body_Y = [l_body(1); 0; l_body(3)];
l_body_Z = [l_body(1); l_body(2); 0];
d_body = B2E\d;
d_body_Y = [d_body(1); 0; d_body(3)];
d_body_Z = [d_body(1); d_body(2); 0];
u_body_Y = [0; 0; UY];
u_body_Z = [0; UZ; 0];
    
% Y plane
X_dot = 0;
X = X + X_dot*ts;
Y_dot = 1/J(2,2) * (cross(l_body_Y,cp) + cross(d_body_Y,cp) + cross(u_body_Y,cu));
Y_dot = Y_dot(2);
Y = (Y + Y_dot*ts);
Z_dot = 1/J(3,3) * (cross(l_body_Z,cp) + cross(d_body_Z,cp) + cross(u_body_Z,cu));
Z_dot = Z_dot(3);
Z = (Z + Z_dot*ts);

% B2E matrices
BX = [1, 0, 0; 0, cos(X), -sin(X); 0, sin(X), cos(X)];
BY = [cos(pi/2+Y), 0, sin(pi/2+Y); 0, 1, 0; -sin(pi/2+Y), 0, cos(pi/2+Y)];
BZ= [cos(Z), -sin(Z), 0; sin(Z), cos(Z), 0; 0, 0, 1];
B2E = BX*BY*BZ; % order of rotations: X, Y, Z

% projection of velocity vector in planes perpendicular to Y and Z body
v_body = B2E\v_unit;
v_body_Y = [v_body(1); 0; v_body(3)];
v_Y = B2E*v_body_Y;
v_body_Z = [v_body(1); v_body(2); 0];
v_Z = B2E*v_body_Z;

% phi angles in planes perpendicular to Y and Z body
phiY_body_dot = (atan(v_Y(1)/v_Y(3)) - phiY_body)/ts;
phiZ_body_dot = (atan(v_Z(2)/v_Z(3)) - phiZ_body)/ts;

xdot(1) = Z_dot;
xdot(2) = phiZ_body_dot;

end