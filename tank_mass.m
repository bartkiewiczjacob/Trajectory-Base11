function [tank] = tank_mass(prop_struct, diameter, material, max_q, max_hydro_press, tank)

%% Function Information 
% Description: 
% Inputs:
%   all structure units in inches and pounds 
%   exception is propellant densities, in lb/ft^3
%   
%   
% Outputs:
% 
% 

%% Constants 
tank(1).name = 'fuel tank';
tank(2).name = 'ox tank';

ullage_perc = 10; % % 
unused_perc = 5; % %
caps = 62; %in^3
pump_LOx_XTRA = 5; % %
pump_LNG_EXTRA = 8; % %

structural_margin = 0.2;
weld_knockdown = .5;

dome_h = diameter/2 / sqrt(2); %square rooted 2 dome height (radial)

tank(1).diameter = diameter;
tank(2).diameter = diameter;
tank(1).dome_h = dome_h;
tank(2).dome_h = dome_h;


%% Calculation
ullage_perc = ullage_perc / 100;
unused_perc = unused_perc / 100;
pump_LOx_XTRA = pump_LOx_XTRA / 100;
pump_LNG_EXTRA = pump_LNG_EXTRA / 100;

%Volume Calculation
volume_mult_fuel = ullage_perc + unused_perc + pump_LNG_EXTRA + 1;
volume_mult_ox = ullage_perc + unused_perc + pump_LOx_XTRA + 1;

tank(1).volume = prop_struct(1).mass / prop_struct(1).density * 12^3 * volume_mult_fuel; %in^3
tank(2).volume = prop_struct(2).mass / prop_struct(2).density * 12^3 * volume_mult_ox; %in^3

tank(1).h = (tank(1).volume - 4/3*pi*diameter^2/4*dome_h)/(pi*diameter^2/4);
tank(2).h = (tank(2).volume - 4/3*pi*diameter^2/4*dome_h)/(pi*diameter^2/4);

%Thickness Calculation
t_press = (max(max(tank.press)) * diameter/2) / (material.minallow * weld_knockdown) * (structural_margin + 1);
p_hydro = max_hydro_press; %* prop_struct(2).density/12^3 * 32.2 * tank(2).h + tank_pressure;

t_hydro = p_hydro * diameter/2 / (2*material.minallow * weld_knockdown)* (structural_margin + 1);

t_required = max(max(t_hydro, t_press));

tank(1).thick = t_required;
tank(2).thick = t_required;

%Cylinder buckling check
f_axial = max_q * pi()*diameter^2/4;

phi = 1/16 * sqrt(diameter/2/t_required);
gam = 1 - .901*(1-exp(-phi));

f_axial_cr = 2*pi*material.E*t_required^2* (gam / sqrt(3-3*material.mu^2)) + min(min(tank.press)) * pi *diameter^2/4;
f_axial_cr = f_axial_cr * .6; %eigenvalue buckling knockdown

if(f_axial_cr <= f_axial)
    error('axial load exceeds critical buckling load')
else
    tank(1).buckling = 'passes cyl buckling check';
    tank(2).buckling = 'passes cyl buckling check';
end

%Mass Calculation
fuel_cyl_vol = (pi/4*diameter^2 - pi*(diameter/2 - t_required)^2)*tank(1).h;
fuel_dome_vol = 4/3*pi*(diameter/2)^2*(dome_h)-4/3*pi*(diameter/2-t_required)^2*(dome_h-t_required);

tank(1).mat_v = (fuel_cyl_vol + fuel_dome_vol + caps*2);
tank(1).mass = tank(1).mat_v * material.density;

lox_cyl_vol = (pi/4*diameter^2 - pi*(diameter/2 - t_required)^2)*tank(2).h;
lox_dome_vol = 4/3*pi*(diameter/2)^2*(dome_h)-4/3*pi*(diameter/2-t_required)^2*(dome_h-t_required);

tank(2).mat_v = (lox_cyl_vol + lox_dome_vol + caps*2);
tank(2).mass = (tank(2).mat_v * material.density);

tank_mass = (tank(2).mass + tank(1).mass);

end