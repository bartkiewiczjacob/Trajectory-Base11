function [propellant, material] = set_tank_param(mass_fuel, mass_ox, name, fuel_density)
% set up for propellant structure
propellant(1).name = 'LNG';
propellant(2).name = 'LOx';

propellant(1).density = fuel_density;
propellant(2).density = 71.24; %lb/ft^3

propellant(1).mass = mass_fuel; %lbm
propellant(2).mass = mass_ox; %lbm

material(1).name = '2195, Sheet';
material(1).density = .098;
material(1).yeild = 66000; %psi
material(1).ult = 70000; % psi
material(1).minallow = 43400; %psi
material(1).E = 1e7; %psi
material(1).mu = .3; %psi


end