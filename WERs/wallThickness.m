function [wall] = wallThickness(p, id)

%% Description

% This function returns the thickness of a tank, given the tank pressure
% and its internal diameter.


%% Initialization

maxStress = 14000; % maximum allowable stress
weldFactor = .85; % welded tank
Y = 0.4;

%% Calculation

wall = p * id / (2 * maxStress * weldFactor + Y * p);

if wall > 0.375
    wall = 0.5;
else
    if wall > 0.25
        wall = 0.375;
    else
        if wall > 0.125
        wall = 0.25;
        else
            wall = 0.125;
        end
    end
end

end
