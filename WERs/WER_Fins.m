% FUNCTION DESCRIPTION
% This function is a weight estimating relationship (WER), for the
% fins on BZB.
%
% ASSUMPTIONS
% For the purposes of generating this weight estimating relationship,
% it was assumed that the fins would have the same geometric shape
% as the fins on BZ. These dimensions are:
%   height:      b_BZ    = 6.16 in
%   tip chord:   c_t_BZ  = 5 in
%   root chord:  c_r_BZ  = 20 in
%   sweep angle: beta_BZ = 66.7 degrees
%   fin area:    S_BZ    = 77 in^2
% The material is assumed to be aluminum with the following
% properties:
%   Shear modulus: G      = 27*10^9 Pa = 3.916*10^6 psi
%   Al Density:    rho_Al = 0.1 lbm/in^3
% For the purposes of fin flutter calculations, the maximum Mach
% number is assumed to be 2.5 at an altitude of 10,000 ft,
% corresponding to the following:
%   Maximum Mach number: M_f = 3
%   Atmospheric pressure: p_atm = 69694.6 Pa
% The following dimensions of BZ are used as a baseline for scaling
% the size of the fins to new rocket sizes:
%   Overall length: L_BZ = 154 in
%   Outer diameter: D_BZ = 6.425 in
% The number of fins is assumed to be three for this weight
% estimating relationship. It may prove advantageous however to use
% four fins.

% INPUTS
% Overall rocket length (inches): L_in
% Outer diameter (inches):        D_in

% OUTPUTS
% Total mass of fins (pound mass): mass_fins_lbm

function [mass_fins_lbm, dimensions] = WER_Fins(L_in, D_in)

    % BZ fin parameters (trapezoid)
    b_BZ_in = 6.16; % in (half-span or fin height)
    c_t_BZ_in = 5; % in (tip chord)
    c_r_BZ_in = 20; % in (root chord)
    S_BZ_in2 = 77; % in^2 (area of fin) 
    % sweep angle = 66.7 deg
    
    % BZ overall parameters
    L_BZ_in = 154; % in (overall length)
    D_BZ_in = 6.425; % in (outer diameter)    
    num_fins = 3;
    
    % Material properties
    G_Al_Pa = 27*10^9; % Pa
    rho_Al_lbmin3 = 0.1; % lbm/in^3
    
    % Flight parameters
    M_f = 3; % Mach number where flutter may occur
    p_atm_Pa = 69694.6; % Pa (Pressure at altitude of max Mach)
    
    % Compute aspect ratio of fin
    AR = (b_BZ_in ^ 2) / S_BZ_in2;
    
    % Compute taper ratio of fin
    lambda = c_t_BZ_in / c_r_BZ_in;
    
    % Compute geometric scaling factor for the fins
    L_bias = 0.75;
    D_bias = 1 - L_bias;
    scale_factor = L_bias*(L_in / L_BZ_in) + D_bias*(D_in / D_BZ_in);
    
    % Compute new root chord, tip chord, and half-span
    c_r_in = c_r_BZ_in * scale_factor; % in
    c_t_in = c_t_BZ_in * scale_factor; % in
    b_in = b_BZ_in * scale_factor; % in
    
    % Compute new fin area
    S_in2 = (b_in / 2) * (c_r_in + c_t_in); % in^2
    
    % Compute thickness required using equation from BZ CDR slide 43
    % Note: this equation likely deserves further investigation
    t_in = ((1.337 * M_f^2 * p_atm_Pa * (1+lambda)) ...
           /(2 * G_Al_Pa * (AR+2)))^(1/3) * (AR * c_r_in);
    
    % Compute mass of fins
    mass_fin_lbm = S_in2 * t_in * rho_Al_lbmin3; % lbm
    mass_fins_lbm = num_fins * mass_fin_lbm; % lbm
    
    % Create vector of fin dimensions for output
    dimensions = [c_r_in, c_t_in, b_in, t_in];
    
end