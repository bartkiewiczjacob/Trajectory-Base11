%% Description

% This executive function estimates the weight of common bulkhead, separate
% bulkhead, and Co-Ax tanks and plots the total weight of the tanks with
% respect to the diameter of tanks.


%% Initialization

% Inputs tank pressure, psi
prompt = '\nEnter the tank pressure (psi): \n';
p = input(prompt); % pressure in psi
while p < 200 || p > 1000
    fprintf('Invalid input: Please try again');
    prompt = '\nEnter the tank pressure (psi): \n';
    p = input(prompt); % pressure in psi
end

mLOX = 25.25; % required mass of LOX, lbm
mCH4 = 8.7; % required mass of methane, lbm
wProp = mLOX + mCH4; % total propellant weight, lbm

prompt = '\nEnter the inner diameter (in) (for CoAx, the ID of LOX tank): \n';
id = input(prompt); % tank ID (for coax, LOX tank ID), in


%% Calculations

fprintf("Total Propellant Weight: %.2f lb\n", wProp);

wCom = calcCommonBulkhead(id, p, mCH4, mLOX);
fprintf("\nCommon Bulkhead Tank Weight: %.2f lb\n", wCom);
fprintf("Common Bulkhead Tank and Propellant Weight: %.2f lb\n", wCom + wProp);

wSep = calcSeparateBulkhead(id, p, mCH4, mLOX);
fprintf("\nCommon Bulkhead Tank Weight: %.2f lb\n", wSep);
fprintf("Common Bulkhead Tank and Propellant Weight: %.2f lb\n", wSep + wProp);

wCoAx = calcCoAx(id, p, mCH4, mLOX);
fprintf("\nCo-Ax Tank Weight: %.2f lb\n", wCoAx);
fprintf("Co-Ax Tank and Propellant Weight: %.2f lb\n", wCoAx + wProp);


%% Plots

% Internal Diameter (in) with respect to Tank Weight (lb) at input pressure
% (psi)
arrayID = 6:.01:8;
arrayIDCom = calcCommonBulkhead(arrayID, p, mCH4, mLOX);
arrayIDTotCom = arrayIDCom + wProp;
arrayIDSep = calcSeparateBulkhead(arrayID, p, mCH4, mLOX);
arrayIDTotSep = arrayIDSep + wProp;
arrayIDCoAx = calcCoAx(arrayID, p, mCH4, mLOX);
arrayIDTotCoAx = arrayIDCoAx + wProp;

figure(1)
subplot(2,1,1)
hold on
plot(arrayID, arrayIDCom, 'b')
plot(arrayID, arrayIDSep, 'r')
plot(arrayID, arrayIDCoAx, 'k')
title(['Inner Diameter (in) vs. Tank Weight (lb) at at ', num2str(p), ' psi'])
xlabel('Inner Diameter (in)')
ylabel('Weight (lb)')
legend('Common Bulkhead', 'Separate Bulkhead', 'Co-Ax', 'location', 'best')
grid on

subplot(2,1,2)
hold on
plot(arrayID, arrayIDTotCom, 'b')
plot(arrayID, arrayIDTotSep, 'r')
plot(arrayID, arrayIDTotCoAx, 'k')
title(['Inner Diameter (in) vs. Tank and Propellant Weight (lb) at ', num2str(p), ' psi'])
xlabel('Inner Diameter (in)')
ylabel('Weight (lb)')
legend('Common Bulkhead', 'Separate Bulkhead', 'Co-Ax', 'location', 'best')
grid on

% Tank Pressure (psi) with respect to Tank Weight (lb) at input internal
% diameter (in)
arrayP = 200:10:1000;
arrayPCom = calcCommonBulkhead(id, arrayP, mCH4, mLOX);
arrayPTotCom = arrayPCom + wProp;
arrayPSep = calcSeparateBulkhead(id, arrayP, mCH4, mLOX);
arrayPTotSep = arrayPSep + wProp;
arrayPCoAx = calcCoAx(id, arrayP, mCH4, mLOX);
arrayPTotCoAx = arrayPCoAx + wProp;

figure(2)
subplot(2,1,1)
hold on
plot(arrayP, arrayPCom, 'b')
plot(arrayP, arrayPSep, 'r')
plot(arrayP, arrayPCoAx, 'k')
title(['Tank Pressure (psi) vs. Tank Weight (lb) at ID = ', num2str(id), ' in'])
xlabel('Tank Pressure (psi)')
ylabel('Weight (lb)')
legend('Common Bulkhead', 'Separate Bulkhead', 'Co-Ax', 'location', 'best')
grid on

subplot(2,1,2)
hold on
plot(arrayP, arrayPTotCom, 'b')
plot(arrayP, arrayPTotSep, 'r')
plot(arrayP, arrayPTotCoAx, 'k')
title(['Tank Pressure (psi) vs. Tank and Propellant Weight (lb) at ID = ', num2str(id), ' in'])
xlabel('Tank Pressure (psi)')
ylabel('Weight (lb)')
legend('Common Bulkhead', 'Separate Bulkhead', 'Co-Ax', 'location', 'best')
grid on