function [total_mass, power_total] = PumpSizingModel(D)
%%
%
% Call: [mass,T,OD,loss] = PumpModel(rpm,D);
%
% Input:
%
%  D       = design specs parameters structure
%   D.Pc    = chamber pressure [psi]
%   D.ox = oxidizer/fluid (string, methane, oxygen, rp-1, or ethanol)
%   D.fuel = fuel/fluid (string, methane, oxygen, rp-1, or ethanol)
%   D.m_dot_ox = flow rate [lbm/s]
%   D.m_dot_fuel = flow rate [lbm/s]
%   D.tb = burn time [s]
%
%  mass     = pump mass [kg]
%  T       = desired torque [N.m]
%  OD       = outer diameter [m]
%  loss     = pump losses [W]
rpm = D.rpm;
assumed_eff = 0.5;
flow_coeff = 0.1;
head_coeff = 0.6;
P0 = 50; %Tank pressure in psi
%Densities assume saturation at 50 psi for cryos and T = 295 K for
%non-cryos
if(D.ox == "oxygen")
    rho_ox = 66.82990;
    P2_ox = D.Pc /0.6;
end
if(D.fuel == "methane")
    rho_fuel = 24.70619;
    P2_fuel = D.Pc /0.6 + 200;
elseif(D.fuel == "rp-1")
    rho_fuel = 50.56665;
    P2_fuel = D.Pc /0.6 + 200;
elseif(D.fuel == "ethanol")
    rho_fuel = 49.19672;
    P2_fuel = D.Pc /0.6 + 200;
end
%%Ox pump. Single stage
iterations = 12;
index = 1;
pump_eff_ox = assumed_eff;
Q_dot_ox_i = D.m_dot_ox / rho_ox; % lbm/s / lbm/ft3 = ft3/s
Q_dot_ox = Q_dot_ox_i * 1.05; % assumes a 2% leak (high estimate)
head_ox = (P2_ox - P0) * 144 / rho_ox; %Rho is in lbm/ft3 so because of the imperial units things we do not include g0
Ns_ox = rpm * ((Q_dot_ox * 448.8) ^ 0.5) / (head_ox ^ 0.75);
pump_eff_list = zeros(iterations,1);
%fprintf("Ns = %f \n", Ns)
while(index <= iterations)
    U2_ox = (head_ox * 32.2 / (pump_eff_ox * head_coeff)) ^ (0.5); %tip speed in ft/s
    D2_ox = U2_ox / (rpm * 2 * pi / 60) * 2; %impeller diameter in ft
    %%FROM NASA SP-8109 Figure 6. Efficiency equations at impeller diameter
    %%inches. Adjusted by 0.08 to match efficiency estimation from mentors for
    %%Ns= 1500 and Diam 2 inches to = 0.5
    c1 = 1.8690476*10^(-10) * D2_ox^4 - 4.47619 * 10 ^(-9) * D2_ox^3 + 3.3792 * 10 ^(-8) * D2_ox^2 - 9.0345238 * 10 ^ (-8) * D2_ox - 4.222857 * 10^-8; - 1.516 * 10 ^(-7);
    c2 = -7.421 * 10 ^(-7) * D2_ox^4 + 1.59353 * 10 ^(-5) * D2_ox^3 - 0.00011386 * D2_ox^2 + 0.000294397 * D2_ox +  0.00016220126;
    c3 = 0.000599458 * D2_ox^4 - 0.012457 * D2_ox^3 + 0.08588 * D2_ox^2 - 0.189068 * D2_ox + 0.396311;
    %%using curve fit above to calc new eff

    pump_eff_ox = c1 * Ns_ox^2 + c2 * Ns_ox + c3; 
   
    power_ox = (rho_ox * head_ox * Q_dot_ox / 550) / pump_eff_ox;
    torque_ox = 63025 * power_ox / rpm;
    %NASA paper with document ID 19750003130, Centrifugal pumps for rocket engines
    %gives a relationship between turbopump weight and the torque of the pump
    %fprintf("\nIteration # %f \n", index)
    %fprintf("Impeller Diameter = %f in, pump_eff = %f, pumpMass = %f lbm\n", D2*12, pump_eff, pumpMass)
    %fprintf("Pump Power = %f hp, Pump Speed = %f rpm, Tip Speed = %f ft/s\n", power, torque, U2)
    pump_eff_list(index) = pump_eff_ox;
    if(index>2)
        if(round(pump_eff_ox,6) == round(pump_eff_list(index - 1),6))
            break
        elseif(round(pump_eff_ox,6) == round(pump_eff_list(index - 2),6))
            %fprintf("Averaging Pump Eff. Old Effs = %f, %f\n", pump_eff, pump_eff_list(index - 1))
            pump_eff_ox = (pump_eff_ox + pump_eff_list(index - 1)) / 2;
            %fprintf("Averaging Pump Eff. New Eff = %f\n", pump_eff)
            U2_ox = (head_ox * 32.2 / (pump_eff_ox * head_coeff)) ^ (0.5); %tip speed in ft/s
            D2_ox = U2_ox / (rpm * 2 * pi / 60) * 2; %impeller diameter in ft
            power_ox = (rho_ox * head_ox * Q_dot_ox / 550) / pump_eff_ox;
            torque_ox = 63025 * power_ox / rpm;
            pumpMass = 5.26 * (torque_ox / 12) ^ 0.638; %mass in lbm
            %fprintf("\nIteration # %f v2 \n", index)
            %fprintf("Impeller Diameter = %f in, pump_eff = %f, pumpMass = %f lbm\n", D2*12, pump_eff, pumpMass)
            %fprintf("Pump Power = %f hp, Pump Speed = %f rpm, Tip Speed = %f ft/s\n", power, torque, U2)
            break
        else
            index = index + 1;
        end
    else
        index = index + 1;
    end
    
end

%%Methane, multistage pump
pump_eff_fuel = assumed_eff;
Q_dot_fuel_i = D.m_dot_fuel / rho_fuel; % lbm/s / lbm/ft3 = ft3/s
Q_dot_fuel = Q_dot_fuel_i * 1.05; % assumes a 2% leak (high estimate)
head_fuel = (P2_fuel - P0) * 144 / rho_fuel; %Rho is in lbm/ft3 so because of the imperial units things we do not include g0
head_fuel1 = head_fuel/2;
head_fuel2 = head_fuel/2;
Ns_fuel1 = rpm * ((Q_dot_fuel * 448.8) ^ 0.5) / (head_fuel1 ^ 0.75);
Ns_fuel2 = rpm * ((Q_dot_fuel * 448.8) ^ 0.5) / (head_fuel2 ^ 0.75);

%fprintf("Ns = %f \n", Ns)
%%Stage1
pump_eff_list = zeros(iterations,1);
iterations = 12;
index = 1;
while(index <= iterations)
    U2_fuel1 = (head_fuel1 * 32.2 / (pump_eff_fuel * head_coeff)) ^ (0.5); %tip speed in ft/s
    D2_fuel1 = U2_fuel1 / (rpm * 2 * pi / 60) * 2; %impeller diameter in ft
    %%FROM NASA SP-8109 Figure 6. Efficiency equations at impeller diameter
    %%inches. Adjusted by 0.08 to match efficiency estimation from mentors for
    %%Ns= 1500 and Diam 2 inches to = 0.5
    c1 = 1.8690476*10^(-10) * D2_ox^4 - 4.47619 * 10 ^(-9) * D2_ox^3 + 3.3792 * 10 ^(-8) * D2_ox^2 - 9.0345238 * 10 ^ (-8) * D2_ox - 4.222857 * 10^-8; - 1.516 * 10 ^(-7);
    c2 = -7.421 * 10 ^(-7) * D2_ox^4 + 1.59353 * 10 ^(-5) * D2_ox^3 - 0.00011386 * D2_ox^2 + 0.000294397 * D2_ox +  0.00016220126;
    c3 = 0.000599458 * D2_ox^4 - 0.012457 * D2_ox^3 + 0.08588 * D2_ox^2 - 0.189068 * D2_ox + 0.396311;
    %%using curve fit above to calc new eff

    pump_eff_fuel = c1 * Ns_fuel1^2 + c2 * Ns_fuel1 + c3; 
   
    power_fuel1 = (rho_fuel * head_fuel1 * Q_dot_fuel / 550) / pump_eff_fuel;
    torque_fuel1 = 63025 * power_fuel1 / rpm;
    %NASA paper with document ID 19750003130, Centrifugal pumps for rocket engines
    %gives a relationship between turbopump weight and the torque of the pump
    %fprintf("\nIteration # %f \n", index)
    %fprintf("Impeller Diameter = %f in, pump_eff = %f, pumpMass = %f lbm\n", D2*12, pump_eff, pumpMass)
    %fprintf("Pump Power = %f hp, Pump Speed = %f rpm, Tip Speed = %f ft/s\n", power, torque, U2)
    pump_eff_list(index) = pump_eff_fuel;
    if(index>2)
        if(round(pump_eff_fuel,6) == round(pump_eff_list(index - 1),6))
            break
        elseif(round(pump_eff_fuel,6) == round(pump_eff_list(index - 2),6))
            %fprintf("Averaging Pump Eff. Old Effs = %f, %f\n", pump_eff, pump_eff_list(index - 1))
            pump_eff_fuel = (pump_eff_fuel + pump_eff_list(index - 1)) / 2;
            %fprintf("Averaging Pump Eff. New Eff = %f\n", pump_eff)
            U2_fuel1 = (head_fuel1 * 32.2 / (pump_eff_fuel * head_coeff)) ^ (0.5); %tip speed in ft/s
            D2_fuel1 = U2_fuel1 / (rpm * 2 * pi / 60) * 2; %impeller diameter in ft
            power_fuel1 = (rho_fuel1 * head_fuel1 * Q_dot_fuel / 550) / pump_eff_fuel;
            torque_fuel1 = 63025 * power_fuel1 / rpm;
            %fprintf("\nIteration # %f v2 \n", index)
            %fprintf("Impeller Diameter = %f in, pump_eff = %f, pumpMass = %f lbm\n", D2*12, pump_eff, pumpMass)
            %fprintf("Pump Power = %f hp, Pump Speed = %f rpm, Tip Speed = %f ft/s\n", power, torque, U2)
            break
        else
            index = index + 1;
        end
    else
        index = index + 1;
    end
    
end

torque_fuel = torque_fuel1 * 2;
torque_total = torque_fuel + torque_ox; 


%%Line and valving mass
A1_ox = Q_dot_ox / 5; %ft3/s / ft/s using a low line velocity, worried about cavitiation
A2_ox = Q_dot_ox / 20; %ft3/s / ft/s using a high line velocity, dont care about cavitation
D1 = 2.875; %outer diameter upstream of the pump
t1 = 0.083; %thickness
D2 = 1.9; %inner diameter downstream of the pump
rho_ss = .289; %lbm/in3 
L1 = 18; %inches
L2 = 18; %inches
if (P2_ox>=1629)
    t2_ox = 0.2;
elseif(P2_ox>=1205)
    t2_ox = 0.145;
elseif(P2_ox>= 705)
    t2_ox = 0.109;
else
    t2_ox = 0.065;
end
if (P2_fuel>=1629)
    t2_fuel = 0.2;
elseif(P2_fuel>=1205)
    t2_fuel = 0.145;
elseif(P2_fuel>= 705)
    t2_fuel = 0.109;
else
    t2_fuel = 0.065;
end

mass_line_ox = ((D1/2)^2 - ((D1/2) - t1)^2) * pi * L1 * rho_ss + ((D2/2)^2 - ((D2/2) - t2_ox)^2) * pi * L2 * rho_ss;
mass_line_fuel = ((D1/2)^2 - ((D1/2) - t1)^2) * pi * (L1) * rho_ss + ((D2/2)^2 - ((D2/2) - t2_fuel)^2) * pi * (L2 +14) * rho_ss;
mass_valve = 10; %lbm

%fprintf("Mass line 1 = %f, Mass line 2 = %f\n", mass_line1, mass_line2)

pumpMass_total = 5.26 * (torque_total / 12) ^ 0.638 + mass_line_ox + mass_line_fuel + 2 * mass_valve; %mass in lbm
mass = (pumpMass_total) / 2.20462; %lbm to kg 
T = torque_total * 0.113; %in*lb to N*m
power_total = power_ox + power_fuel1 * 2;
motor_mass = power_total /9.3;
battery_mass = 12.76 * power_total * 1.1 * (D.tb/60/60);
inverter = 10;
total_mass = pumpMass_total + motor_mass + battery_mass + inverter;

loss = power_ox * 745.7 * pump_eff_ox; %Loss. Power is converted from hp to w then mutliplied by the pump eff
if D2_ox > D2_fuel1
    OD = D2_ox * 1.5;
else
    OD = D2_fuel1 * 1.5;
end
