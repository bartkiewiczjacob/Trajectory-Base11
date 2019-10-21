function [tank, helium, material,max_hydro_press] = calc_tanks(Pc,pressfed,tank_temp,diameter,max_q,mass_fuel,mass_ox, oxidizer_mass, acceleration_g, fuel_density)

tank = tank_press_calc(Pc,tank_temp,pressfed);
max_hydro_press = hydro_accel_calc(acceleration_g,oxidizer_mass,diameter);
[propellant, material] = set_tank_param(mass_fuel, mass_ox, fuel_density);
tank = tank_mass(propellant, diameter, material(1), max_q, max_hydro_press, tank);
[helium, tank] = helium_mass(Pc, tank, pressfed, tank_temp);

end