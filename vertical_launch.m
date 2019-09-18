function [dhdt] = vertical_launch(t, h, thrust, wmass, mdot, burn_time)

if t < burn_time
    T = (thrust)*4.44822; % N
else
    T = 0;
end
m = (wmass - mdot*t)*0.453592; % kg

Cd = 0.3;
d = 18*0.0254; % m
S = pi*d^2/4; % m2
G = m*9.81; % N

[temp, A, P, RHO] = atmosisa(h(1));

dhdt = zeros(2,1);
dhdt(1) = h(2);
dhdt(2) = 1/m*(-Cd*0.5*(RHO)*h(2)*S + T - G);




