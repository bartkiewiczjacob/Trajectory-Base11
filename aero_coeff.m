function [CD, CL, CP] = aero_coeff(Ma, alpha)
% compute the drag and lift coefficients and the center of pressure for
% given Mach number and angle of attack

% access aerodynamics test data
global Ma_val CD_data_mach CL_data_mach CP_data_mach CD_fit_alpha CL_fit_alpha CP_fit_alpha time

% find the coefficients from Ma (alpha = 0)
index = find(abs(Ma_val-Ma) <= 1e-10);
CD_0 = CD_data_mach(index);
CL_0 = CL_data_mach(index);
CP_0 = CP_data_mach(index);

[val, index] = min(abs([0.1, 0.5, 1.1, 2, 5]-Ma)); % check which Ma case is closest
% interpolate from the alpha = 0 value using the polyfit with closest Ma
CD = CD_fit_alpha(1, index)*alpha^2 + CD_fit_alpha(2, index)*alpha + CD_0;
CL = CL_fit_alpha(1, index)*alpha + CL_0;
CP = CP_fit_alpha(1, index)*alpha^4 + CP_fit_alpha(2, index)*alpha^3 + ...
    CP_fit_alpha(3, index)*alpha^2 + CP_fit_alpha(4, index)*alpha + CP_0;

% if time == 25.5
%     keyboard
% end
end
