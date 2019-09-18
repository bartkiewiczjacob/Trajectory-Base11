clear; clc;
T = 5000*4.44822;
Cd = 0.3;
d = 18*0.0254; % m
S = pi*d^2/4; % m2
burn_time = 40;
wmass = 1250*0.453592;
mdot = 18.75*0.453592;

sim('test.slx')


%%
clear; clc;
thrust = 5000;
wmass = 1250;
mdot = 18.75;
burn_time = 40;
[t,h] = ode45(@(t,h) vertical_launch(t, h, thrust, wmass, mdot, burn_time), [0, 300], [0 0 wmass*0.453592])
max(h(:,1))