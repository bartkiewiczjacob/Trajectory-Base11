function [wTot, odCH4, tCH4, hTot, hLOX, hCH4] = tankWER(id, p, Isp, OF)

%% Description

% This function estimates the weight of common bulkhead, separate
% bulkhead, and Co-Ax tanks and plots the total weight of the tanks with
% respect to the diameter of tanks.


%% Initialization

fCooling = 1.1; % film cooling percentage
% mLOX = 32.6676; % required mass of LOX, lbm
% mCH4 = 9.22 * fCooling; % required mass of methane, lbm

impulse = 9208; % lbf*s

mTot = impulse / Isp; % lb

mLOX = (mTot/(OF+1))*OF;
mCH4 = (mTot/(OF+1));
mFilm = mCH4*(fCooling-1);
mCH4 = mCH4*fCooling;

%% Calculations

% [wTot, odCH4, odLOX, tCH4, tLOX, hCH4, hLOX] = calcCoAx(id, p, mCH4, mLOX);
[wTot, odCH4, tCH4, hTot, hLOX, hCH4] = calcCommonBulkhead(id, p, mCH4, mLOX);
wTot = wTot + mFilm;
