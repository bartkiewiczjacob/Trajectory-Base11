function [helium,tank] = helium_mass(Pc, tank, pressfed, tank_temp)

if pressfed == 1
    
%     tank_temp = -150; %F
%     tank(1).press  = Pc / .6; %kyle approx. fuel tank
%     tank(2).press = tank(1).press + 200; %lox tank + regen
     R = 40.8829;

    fuel_he_moles = tank(1).volume*tank(1).press / (R*(tank_temp+359.67));
    lox_he_moles = tank(2).volume*tank(2).press / (R*(tank_temp+359.67));

    he_moles = fuel_he_moles + lox_he_moles;
    helium.mass = he_moles/.24983748071879*.00220462;

    helium.temp = 70; %F
    helium.press = 4500; %psi
    helium.gas_vol = he_moles * R*(helium.temp + 459.67)/helium.press;

    helium.tank = helium.gas_vol * helium.press / 1e6; %lbm

    helium.total_mass = helium.tank + helium.mass;
elseif pressfed == 0
    helium.total_mass = 20; %lbm
end

end

