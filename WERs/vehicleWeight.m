function [mass] = vehicleWeight(diameter, height)

    % height and diameter in inches
    % calculates hollow cylinder mass

    density = .1;    %lbs/in^3
    wallThickness = .125; %in
    volume = pi() * ((diameter / 2)^2 - (diameter/2 - wallThickness)^2) * height;
    mass = volume * density;
end