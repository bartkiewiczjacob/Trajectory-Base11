clc

maxVelocity = 2250; %ft/s
altitude = 10000; %ft

[rho,speedOfSound,T,pressure,nu,z] = atmos(altitude, 'units', 'US');

% speed of sound in ft/s, pressure in lb/ft^2

pressure = pressure / 144; %convert to psi
seaPressure = 14.69594878; %psi

tipChord = 5.44;  %in
rootChord = 24; %in
avgChord = (tipChord + rootChord) / 2;
semispan = 17.5; %in 

G = 3770000; %psi
AR = (semispan^2) / (0.5 * (rootChord + tipChord) * semispan); %unitless
lambda = tipChord / rootChord; %unitless

for thickness = 0.1:0.01:5 

    term2 = ((39.3*(AR^3)/(((thickness/avgChord)^3)*(AR+2))) * ((lambda+1)/2) * (pressure/seaPressure));
    flutter2 = speedOfSound * sqrt(G/term2);
       
    if (flutter2 > 4343 * 1.5)
        thickness
        flutter2
    end
        
end

