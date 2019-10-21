function [sol] = Bending(qload, L, b, h, type)

%% Initial calualtions
V = qload*L;        % Shear force      
M = qload*L/2;     % Max moment

if type == 1
    Ic = b*h^3/12;   % Centroidal Moment
else
    Ictest = b*h^3/12;
    Af = 0.5*h;
    Ic = (0.0432/12*h^3 + 2*(0.0048/12*h^3 + (0.3*h)^2*Af))*type;
end
Q = b*h^2/8;        % First Moment of Area 

A = h*b;            % Largest area for shear stress

%% Stress, Strain, and Deflection

sigma = M*h/2/Ic;
tau = 3*V/2/A;

sol = [sigma, tau];



