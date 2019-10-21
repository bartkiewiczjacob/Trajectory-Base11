function [T_h] = getEngineAtAlt(T_opt,AR,P_h,Pc,At,Pe)
%% Engine vs Altitude
%
%   Inputs:
%       T_opt       (Optimum thrust)                [Lb-f]
%       isp_opt     (Optimum Isp                    [s]
%       AR          (Area ratio)                    [-]
%       cstar       (Characteristic Velocity)       [ft/s]
%       P_h         (Pressure at altitude h)        [Psi]
%       At          (Throat Area)                   [in2]
%       Pe          (Exit Pressure)                 [psi]
%
%   Outputs: 
%       T_h         (Thrust at given altitude)      [Lb-f]
%       isp_h       (Isp at a given altitiude)      [s]
%

%% Constants
g = 32.2; % [ft/s^2]

%% Calculations
T_h = T_opt + (Pc.*At.*AR.*((Pe-P_h)/Pc));
end
