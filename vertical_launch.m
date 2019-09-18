function [dhdt] = vertical_launch(t, h, thrust, wmass, mdot, burn_time)

if t < burn_time
    T = (thrust)*4.44822; % N
    mdot = mdot*0.453592;
else
    T = 0;
    mdot = 0;
end

Cd = 0.3;
d = 18*0.0254; % m
S = pi*d^2/4; % m2
G = h(3)*9.81; % N

if h(1) < 84000
    [temp, A, P, RHO] = atmoscoesa(h(1),'None');
else
    RHO = 0;
end

dhdt = zeros(2,1);
dhdt(1) = h(2);
dhdt(2) = 1/h(3)*(-Cd*0.5*(RHO)*h(2)^2*S + T - G);
dhdt(3) = -mdot;


