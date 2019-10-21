function [propellant, material] = prop_params()
% set up for propellant structure

propellant(1).name = 'LNG';
propellant(2).name = 'LOx';

propellant(1).density = 26.37; %lb/ft^3
propellant(2).density = 71.24; %lb/ft^3

propellant(1).mass = 217.03; %lbm
propellant(2).mass = 585.97; %lbm

material(1).name = '2195, Sheet';
material(1).density = .098;
material(1).yeild = 66000; %psi
material(1).ult = 70000; % psi
material(1).minallow = 43400; %psi
material(1).E = 1e7; %psi
material(1).mu = .3; %psi


end