%Assumptions: the rocket cannot weigh more than 100-something lbs because
%of physics. As a result, tube thickness is governed more by manufacturing
%and availability concerns than structural ones. It will be difficult to
%work with tubing much more than .125", and that thickness will almost
%certainly be enough to sustain any reasonable loading with the size of
%rocket we're expecting. 

% Takes vehicle diameter and length and returns the actual diameter used
% (limited by available stock) and the weight of the shell. 
% Material is 0 for aluminum or 1 for fiberglass
% All units in inches and lbs
function [diameter, weight] = Vehicle_WER(D, L, material)
    existing = linspace(4,10,25);   %Readily available sizes from https://www.vitaneedle.com/aluminum-tube/

    density = .1;    %lbs/in^3
    wallThickness = .125; %in
    
    if material == 1
        existing = [4.5 5.15 5.375 6.17 7.518 8.005];
        ID = [4.375 5 5.525 6 7.708 7.815];
        linearDensity = [16.8 19 16.8 24 25 27]./16./12; %Converting from oz/ft to lb/in
    end
    [M,I] = min(abs(existing - D));
    diameter = existing(I);
    if material == 1
%         wallThickness = (diameter - ID(I))/2;
        weight = linearDensity(I) * L;
        return;
    end
    
    volume = pi() * ((diameter / 2)^2 - (diameter/2 - wallThickness)^2) * L;
    weight = volume * density;
end