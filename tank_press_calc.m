function [tank] = tank_press_calc(Pc,tank_temp,pressfed)

if pressfed == 1
    tank_temp = -150; %F
    tank(1).press  = Pc / .8; %kyle approx. fuel tank
    tank(2).press = tank(1).press + 200; %lox tank + regen
else
    tank(1).press = 50;
    tank(2).press = tank(1).press;
    
end
