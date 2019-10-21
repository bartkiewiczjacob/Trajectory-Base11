format compact; clc; clear all; close all;

max_q = 17.93;
tank_presssure = 100
diameter = 18;
acceleration_array = [1.5,1.65,1.96];
oxidizer_m_array = [970, 563,331.78];

mass_fuel = 800/2.7;
mass_ox = 800 - 800/2.7;

[propellant, material] = set_tank_param(mass_fuel, mass_ox);

max_hydro_p = hydro_accel_calc(acceleration_array, oxidizer_m_array, diameter);
tank = tank_mass(propellant, tank_presssure, diameter, material(1), max_q,max_hydro_p)
intertank = intertank_mass(tank,material(1))

