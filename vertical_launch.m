function [dhdt] = vertical_launch(t, h, thrust, mdot, burn_time, Cd, De, Dt, Pc, Pe, S)
if h(1) < 84000
    [temp, A, P, RHO] = atmoscoesa(h(1),'None');
else
    RHO = 0;
    P = 0;
end
if t < burn_time
    T = getEngineAtAlt(thrust, De^2/Dt^2, P*0.000145038, Pc, pi*Dt^2/4, Pe);
    T = T*4.45; %N
    mdot = mdot*0.453592;
else
    T = 0;
    mdot = 0;
end

d = 18*0.0254; % m
G = h(3)*9.81; % N
S = pi*d^2/4;
dhdt = zeros(2,1);
dhdt(1) = h(2);
dhdt(2) = 1/h(3)*(-Cd*0.5*(RHO)*h(2)^2*S + T - G);
dhdt(3) = -mdot;



