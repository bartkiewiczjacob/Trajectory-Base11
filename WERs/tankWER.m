function [wTot, odLOX, odCH4, tLOX, tCH4, hLOX] = tankWER(id, p)

%% Description

% This function estimates the weight of common bulkhead, separate
% bulkhead, and Co-Ax tanks and plots the total weight of the tanks with
% respect to the diameter of tanks.


%% Initialization

fCooling = 1.1; % film cooling percentage
mLOX = 32.6676; % required mass of LOX, lbm
mCH4 = 9.22 * fCooling; % required mass of methane, lbm


%% Calculations

[wTot, odCH4, odLOX, tCH4, tLOX, hCH4, hLOX] = calcCoAx(id, p, mCH4, mLOX);
