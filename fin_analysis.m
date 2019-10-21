function [data] = fin_analysis(Diameter, fin_num, mat_prop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% Calculates the fin size and the fin mass. 
% Assuming flat plate theory and honeycomb structure.
% 
% Input: (Diameter, fin_num)
% Diamter => Diamter of rocket in inches
% fin_num => number of fins (4)
%
% Output: data it is made of [mass a]
% mass => the mass of each fin in lbs
% a => root chord of fin in inches
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input
% Diameter = 18;  % Diameter in inches
% fin_num = 4;    % Number of fins

%% Constants and Data 
% Orientation

% Aerodynamics
q_max = 17.93;      % Max dynamic pressure [psi]

%% Dimensions

dims = hist_sizing(Diameter,fin_num);

a = dims(1);             % Long length of Trapazoid [in]
b = dims(2);             % Short length of Trapazoid [in]
h = dims(3);             % Height of Trapazoid  [in]

span_avg = (h*b + h/3*(a-b))/a;
span_max = h;
cord_avg = b+a/2;
cord_max = a;
area = (a+b)/2 *h;
beam_num = floor(cord_avg/2.5);

% Material Properties
sf = 1.2;           % Safety factor
% mat_name = xlsread('Material_prop.xlsx')


%% Drag forces
n = 1;
start = 1;
last =  5;
angles = linspace(start,last, 5);

for alpha = angles
th = 1;
L = q_max*pi*sin(alpha)*cord_avg;
D = q_max*pi*cos(alpha)*cord_avg;
Fd_bend1 = L*cos(alpha) + D*sin(alpha); % Drag force 1
Fd_bend2 = 1/2*q_max*th*cd;             % Drag force 2

%% Thickness viration
thick = 0.1:0.05:0.3;
    
    j = 1;
    warp = [];
    bend = [];
    for th = thick
        sidebend_max = Bending(Fd_bend1, cord_avg, span_avg, th,1);
        topbend_max =  Bending(Fd_bend2, span_max, cord_max, th,1);
        warp(j,:) = topbend_max;
        bend(j,:) = sidebend_max;
        j = j+1;
    end
    
    warp2(:,:,n) = warp;
    bend2(:,:,n) = bend;
    
    n = n+1;
end

%% Thickness
stress_cap = (0.5*mat_prop{2,2} +  0.5*mat_prop{3,2})/sf;
diff = abs(stress_cap - abs(bend2(:,1,4)));
ind = find(diff == min(diff));
ideal_th = ceil(thick(ind)*100)/100;

volume = ideal_th*area;
mass = (.5*mat_prop{2,5} +  .5*mat_prop{3,5})*volume; 

data = [mass a];
%% Output
% fprintf("\nThe mass of the fin would be %.2f lbs with thickness of %.2f"". \n\n",mass, ideal_th)


%% Plot 
% 
% figure
% subplot(2,1,1)
% hold on
% for m = 1:n-1
%     plot (thick, abs(bend2(:,1,m)))
% end
% 
% for m = 1:height(mat_prop)
%     plot([min(thick), max(thick)], [mat_prop{m,2}, mat_prop{m,2}]/sf,'--')
% end
% grid on
% grid minor
% title('\sigma "Normal Stress" z axis')
% xlabel('Thickness [inches]')
% ylabel('\sigma Normal Stress [psi]')
% legend(num2str(angles'))
% 
% subplot(2,1,2)
% hold on
% for m = 1:n-1
%     plot (thick, abs(bend2(:,2,m)))
% end
% grid on
% grid minor
% title ('\tau "Shear stress"')
% xlabel('Thickness [inches]')
% ylabel('\tau Shear stress [psi]')
% legend(num2str(angles'))
% 
% 
% figure
% subplot(2,1,1)
% hold on
% for m = 1:n-1
%     plot (thick, abs(warp2(:,1,m)))
% end 
% 
% for m = 1:height(mat_prop)
%     plot([min(thick), max(thick)], [mat_prop{m,2}, mat_prop{m,2}]/sf,'--')
% end
% grid on
% grid minor
% title('\sigma "Normalstress" x-axis')
% xlabel('Thickness [inches]')
% ylabel('\sigma Normal Stress [psi]')
% legend(num2str(angles'))
% 
% subplot(2,1,2)
% hold on
% for m = 1:n-1
%     plot (thick, warp2(:,2,m))
% end
% grid on
% grid minor
% title ('\tau "Shear stress"')
% xlabel('Thickness [inches]')
% ylabel('\tau Shear stress [psi]')
% legend(num2str(angles'))

%% End it 



