function [weight] = avionicsWER(numAltimeters, numTC, volHousing)
%AVIONICSWEIGHT Calculates estimated weight of Avionics in BZB

% expected values 10-19-19: 
%   numAltimeters = 2
%   numTC = 2
%   volHousing = 30

%% Constants
%  Weight (lb) of components

rpi = 0.0925942;        % Raspberry Pi
mcc118 = 0.0330693;     % MCC 118 DAQ
raven = 0.0145505;      % Featherweight Raven 3 Altimeter
rrc3 = 0.0375;          % RRC3
tcconv = 0.00242508;    % TC Converter
gps = 0.0220462;        % Featherweight GPS
batt = 0.099208;        % 9V Battery
pbatt = 0.17637;        % Portable Charger (Anker Lipstick-size)

pladensity = 0.0289018; % Density of PLA
wire = 0.220462;        % Worst-case weight of wires (100g)

%% Calculations

weightbatteries = (numAltimeters + numTC) * batt;
weightaltimeter = numAltimeters * raven;
weighthousing = volHousing * pladensity;
weighttc = numTC * tcconv;

weight = weightbatteries + weightaltimeter + weighthousing + weighttc ...
    + rpi + mcc118 + rrc3 + gps + pbatt + wire;
end

