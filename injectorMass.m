function injector_Mass = injectorMass(Dc,rho)
%% Injector Mass Estimate
%
%   rho = 0.322978 CuCrZr
Dc_CAD = 5.2;
m_CAD = 6;
Ac_CAD = (0.25*pi*Dc_CAD^2)*(rho/0.0975);
M_SA = m_CAD/Ac_CAD;

Ac = 0.25.*pi.*Dc.^2;

injector_Mass = M_SA.*Ac;

end


