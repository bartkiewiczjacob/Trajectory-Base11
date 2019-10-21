function [data] = Boattail(diamter, L1, d_ne,L2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boat-tail and fin holder cylinder mass
%
% Input: (diamter, L1, r_ne,L2)
% diameter => diamter of the rocket in inches
% L1 => distance between end of tank to end of engine
% r_ne => nozzle exit diamter inchea
% L2 => rootchord of the fins
% 
% Output: data is made of [Mass_b Mass_c]
% Mass_b => mass of the boat-tail in lbs
% Mass_c => mass of cylinder to hold fins in lbs
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Boattail sizing
% clear
% close all
% clc

%% Input 
% L1 = 36;        % Distance from end of tank to end of nozzle [in]
% 
r_ne = d_ne/2;   % Nozzle exit radius [in]
r_o = diamter/2;     % Rocket radius [in]
th = .25/2;       % Thicness of boattail [in] (Assume same as fin)
% L2 = 28.8;      % Root chord of fin [in]                

%% Calculations
Lt = L1+1;  
r_be = r_ne+2.5;  % Exit radius of boattail [in]
h1 = Lt - L2;   % Length of boat tail [in]

ht = Lt*r_o/(r_o-r_be); % Height of boattail if it was a closed cone [in]

Vol_b = pi*(h1)/3 * (2*r_o*th - th^2);
Mass_b = Vol_b*0.1;

Vol_c = pi*L2*(2*r_o*th - th^2);
Mass_c = Vol_c*.1;

data = [Mass_b Mass_c];

%% Output 
% fprintf('At a thickness of %.2fin the boatail would have a length of %.2fin, exit radius %.2fin, and weigh %.2f lbs\n\n',[th,h1,r_be,Mass_b])
% fprintf('At a thickness of %.2fin the cylinder to support fins would have a length of %.2fin and weigh %.2f lbs\n\n',[th,L2,Mass_c])
