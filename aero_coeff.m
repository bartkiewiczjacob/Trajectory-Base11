function [CD, CL, CP] = aero_coeff(Ma, alpha)
global Ma_val CD_data_mach CL_data_mach CP_data_mach CD_fit_alpha CL_fit_alpha CP_fit_alpha

index = find(Ma == Ma_val);
CD_0 = CD_data_mach(index);
CL_0 = CL_data_mach(index);
CP_0 = CP_data_mach(index);

[val, index] = min(abs([0.1, 0.5, 1.1, 2, 5]-alpha));
CD = CD_fit_alpha(1, index)*alpha^2 + CD_fit_alpha(2, index)*alpha + CD_0;
CL = CL_fit_alpha(1, index)*alpha + CL_0;
CP = CP_fit_alpha(1, index)*alpha^4 + CP_fit_alpha(2, index)*alpha^3 + CP_fit_alpha(3, index)*alpha^2 + CP_fit_alpha(4, index)*alpha + CP_0;
end
