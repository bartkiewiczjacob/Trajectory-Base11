function cdb = basedrag(reynolds, wetted_surf, cross_section, l_over_d)
%% BASE DRAG
%
% Author: Ben Kuras
%
% Description: This function takes in rocket geometry and reynolds number
% over the flight and calculates the base drag coefficient
%
% Inputs
% - reynolds: reynolds number
% - wetted_surf: wetted surface area of the rocket (must be same units as
% cross_section)
% - cross_section: cross section of the rocket (must be same units as
% wetted_surf)
%
% outputs
% - cdb: Base drag coefficient
%% INITIALIZATION
cf = 0.027 ./ (reynolds).^(1/7); %Calculates skin friction coefficient assuming turbulent boundary layer
cdn_cdbt = 1.02*cf .* (1+1.5/l_over_d^(3/2)) * wetted_surf/cross_section; %Caclulates sum of nose cone and body tube drag coefficient
cdb = 0.029./sqrt(cdn_cdbt); %Calculates base drag coefficient
