function [mass,isp,isp_vac,cf,cf_vac,Dt,Dc,De,mdot] = engineMass(OF, Pc, Pe, L_star, CR, rho, eta_cstar, eta_cf, thrust, prop, ffc)
%%  Engine Sizing and Mass Estimate Function
%
%   Estimates engine mass based on 
%
%   Inputs:     OF (Mixture ratio) [-]
%               Pc (Chamber Pressure [PSIA]
%               Pe (Design Exit Pressure) [PSIA]
%               L_star (Characteristic Length) [IN]
%               CR (Contraction Ratio) [-]
%               rho (Material density) [LB-M/IN^3]
%               etacstar (combustion efficiency) 
%               etacf (nozzle efficiency)
%               thrust (thrust) [lbf]       
%               Twg (Wall side gas temperature) [R]
%               prop (Propellant Type) KEY: Methane N=1, RP-1 N=2, ETHANOL N=3
%               ffc (% Of fuel that's used for film cooling) [%]
%
%   Outputs:    mass (Engine Mass) [LB-M]
%               isp (Specific Impulse) [S]
%              
%


%% Useful Constants

g = 32.174; % Gravity [ft/s^2]
%% CEA

% Fluid Properties
T_LOX = 90.17; % Liquid Oxygen temperature [K]

if prop == 1
    
    % Methane Fuel CEA Call
    T_LCH4 = 111;   % Liquid Methane Temperature [K]
    engine = CEA('problem','rocket','equilibrium','o/f',OF,'pi/p',Pc/Pe,'case','AAE 450 TCA METHANE','p(psi)',Pc,'reactants','fuel','CH4(L)','wt%',100,'t(k)', T_LCH4,'oxid','O2(L)','wt%',100,'t(k)',T_LOX,'output','eng','end');

elseif (prop == 2)
        
        % RP-1 Fuel CEA Call
        T_RP1 = 300;    % RP-1 Temperature [K]        
        engine = CEA('problem','rocket','equilibrium','o/f',OF,'pi/p',Pc/Pe,'case','AAE 450 TCA RP1','p(psi)',Pc,'reactants','fuel','RP-1','wt%',100,'t(k)', T_RP1,'oxid','O2(L)','wt%',100,'t(k)',T_LOX,'output','eng','end');

    elseif (prop == 3)
            
            % Ethanol
            T_eth = 300;                    % Ethanol Temperature [K]             
            engine =    CEA('problem','rocket','equilibrium','o/f',OF,'case','AAE 450 TCA ETHANOL','p(psi)',Pc, 'pi/p', Pc/Pe, 'reactants','fuel','C2H5OH(L)','wt%',75,'T(k)',T_eth,'fuel','H2O(L)','wt%',25,'T(k)',T_eth,'oxid','O2(L)','wt%',100,'t,k',T_LOX,'output','eng','end');
            
        else
            error('\nYOU MUST CHOOSE A VALID VALUE FOR N!!!\n\nMETHANE: N=1\nRP-1: N=2\nETHANOL: N=3\n')
end
            
% Pull Engine Characteristics from CEA

isp = eta_cstar*eta_cf*engine.output.eql.isp(3);            % Real Specific Impulse at design condition [s]
isp_vac = eta_cstar*eta_cf*engine.output.eql.isp_vac(3);    % Real Vacuum Specific Impulse [s]
cstar = eta_cstar* engine.output.eql.cstar(3);              % Real C* [ft/s]
cf = eta_cf*engine.output.eql.cf(3);                        % Real Thrust Coefficient [-] 
cf_vac = eta_cf*engine.output.eql.cf_vac(3);

AR = engine.output.eql.aeat(3);                             % Expansion Ratio [-]

%% Calculations

mdot = thrust/isp;                 % Mass flow rate [lbm/s]
mdot_f = mdot/(1+OF);               % Fuel Mass Flow [lbm/s]
mdot_f_total = (1+ffc)*mdot_f;      % Total fuel mass flow (Combined with film cooling) [Lbm/s]
mdot_ffc = mdot_f_total - mdot_f;   % Fuel film cooling %
mdot = mdot + mdot_ffc;             % Total mass flow rate adjusted for film cooling

At = cstar*mdot/(Pc*g);              % Throat area [in^2]
Dt = sqrt(4*At/pi);                % Diameter of Throat [in]

Ae = AR*At;                        % Area of Exit [in^2]
De = sqrt(4*Ae/pi);                 % Diameter of Exit [in]

Lc = L_star/CR;                     % Calculates Chamber Length [in]
Ac = CR*At;                         % Chamber Area [in^2]
Dc = sqrt(4*Ac/pi);             % Chamber Diameter [in]
theta_i = 38;
theta_e = 9;

%%  Nozzle & Combustor Mass Estimate ( 80% Bell Nozzle)

% Nozzle Mass Estimate

[x_1,r_1] = getContour(Dt,AR,CR,L_star,theta_i,theta_e);    % Nozzle contour outputs (Axial,Radial) 

SA_Nozzle = 0;
for idx = 1:(length(x_1)-1) 
    dx = abs(x_1(idx+1)-x_1(idx));
    SA_Nozzle = SA_Nozzle + 2.*pi.*dx.*(r_1(idx));
end
Mass_SA_CAD_Model = (46.42/418.91)*(rho/0.321);    % 418.91 [lbm/in^2]
Mass_nozzle = SA_Nozzle*Mass_SA_CAD_Model; % Nozzle Mass [Lbm]

% Injector Mass Estimate
Dc_CAD = 5.2;  
m_CAD = 6;
Ac_CAD = (0.25*pi*Dc_CAD^2)*(rho/0.0975);
M_SA = m_CAD/Ac_CAD;
Ac = 0.25.*pi.*Dc.^2;

injector_Mass = M_SA.*Ac;

mass = Mass_nozzle + injector_Mass+15; %[lbm]
end