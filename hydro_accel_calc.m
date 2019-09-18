function [max_hydro_pressure] = hydro_accel_calc(acceleration_g, oxidizer_mass, diameter)


[max_accel,index] = max(acceleration_g .* oxidizer_mass);

max_load = max_accel;

SA = pi*diameter^2 /4;

max_hydro_pressure = max_load/SA;

end

