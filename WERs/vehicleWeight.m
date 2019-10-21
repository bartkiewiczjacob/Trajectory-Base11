function [weight] = vehicleWeight(diameter, height)
    density = .1;    %lbs/in^3
    wallThickness = .125; %in
    volume = pi() * ((diameter / 2)^2 - (diameter/2 - wallThickness)^2) * height;
    weight = volume * density;
end