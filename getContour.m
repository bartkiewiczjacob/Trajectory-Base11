function [x, r] = getContour(Dt, AR, CR, L_star, initial_parabola_angle, final_parabola_angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code graphs the engine contour assuming 80% bell nozzle - ideal for contour
% (includes combustion chamber, throat, and nozzle)
%
%
%   OUTPUTS:
%           - x (Axial distance) [in]
%           - r (Radial Distance) [in]
%
%   INPUTS: 
%           - Dt (Throat Diameter) [in]
%           - AR (Expansion Ratio) [-] %expansion ratio, used to determine contour of the nozzle design (ratio of the nozzle exit area, to the nozzle throat area)
%           - CR (Contraction Ratio) [-] %ratio of chamber area to throat area
%           - L_star (Chamber Length * CR) [in]
%           - initial_parabola_angle [deg] %approximation, 80% bell nozzle (referencing G.V.R. Rao angle vs. expansion ratio graph)
%           - final_parabola_angle [deg] %approximation, 80% bell nozzle (referencing G.V.R. Rao angle vs. expansion ratio graph)
%
%
% Helpful resources for reference (for Quadratic Bezier Curve geometric explanation): 
% Sutton 8th edition Chapter 3
% https://engineering.purdue.edu/AAECourses/aae439/2008/aae439_class_lecture/lecture_notes/chap4_71_95.pdf
% http://graphics.cs.ucdavis.edu/education/CAGDNotes/Quadratic-Bezier-Curves.pdf
% http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
%
% Written by Caroline Kren on 09/30/2018
% Last Updated by Caroline Kren on 09/17/2019
%
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Notes
%when updating the expansion ratio, update the final and initial parabola angles as well with G.V.R. Rao angle vs. expansion ratio graph pg. 80 8E Sutton
%this code creates contour by creating contour of 1.5 converging throat, 0.4 diverging throat, and parabolic nozzle first, then creating contour for chamber converging and chamber itself

%% Initializations

Rt = Dt/2; %radius at the throat of the nozzle [in]
length_chamber_convergence_section = 3; %length of chamber convergence section

% %TEST CASE USING SUTTON %only use from 1.5 converging throat to end of nozzle
% %constants using for reference (using Sutton example for code case example)
% throat_radius = 1 %radius at the throat of the nozzle 
% contraction_ratio = 8; %ratio of chamber area to throat area
% expansion_ratio = 25; %expansion ratio, used to determine contour of the nozzle design (ratio of the nozzle exit area, to the nozzle throat area)
% initial_parabola_angle = 30; %approximation, 80% bell nozzle (referencing G.V.R. Rao angle vs. expansion ratio graph)
% final_parabola_angle = 8.5; %approximation, 80% bell nozzle (referencing G.V.R. Rao angle vs. expansion ratio graph)

%Sub equations using input parameters
Lc=5;                                       % length of combustion chamber (defined as from injector to start of converging section - straight portion)
throat_area = pi * Rt^2;                                % throat area (inches^2)
chamber_area = CR * throat_area;                        % chamber area
chamber_radius = sqrt(chamber_area / pi);               % chamber radius 
exit_area = AR * throat_area;                           % area at the exit of the nozzle
exit_radius = sqrt(exit_area / pi);                     % radius at the exit of the nozzle
cone_length = (((sqrt(AR)-1)*Rt)/(tand(15)));           % Rao cone length equation
eighty_percent_bell_contour_length = cone_length * 0.8; % 80% bell nozzle length 

%1.5*THROAT_RADIUS CONVERGING PORTION
angle_1_5 = linspace(-135,-90);                         %angles for 1.5 throat radius curve to be graphed
x_1_5 = 1.5 * Rt * cosd(angle_1_5);                     %x as function of angle_1_5
y_1_5 = (1.5 * Rt * sind(angle_1_5)) + (1.5 * Rt) + Rt; %y as function of angle_1_5

%0.4*THROAT_RADIUS DIVERGING PORTION
angle_0_4 = linspace(-90,(initial_parabola_angle-90));  %angles for 0.4 throat radius curve to be graphed
x_0_4 = 0.4 * Rt * cosd(angle_0_4);                     %x as function of angle_0_4
y_0_4 = (0.4 * Rt * sind(angle_0_4)) + (0.4 * Rt) + Rt; %y as function of angle_0_4

%DIVERGING PARABOLA SECTION
%Determining parameters for Quadratic Bezier Curve
%Point P0 (point of inflection between semicircle and parabola)
P0_angle = initial_parabola_angle - 90;                 %angle of P0
P0_x = 0.4 * Rt * cosd(P0_angle);                       %finding x value of P0
P0_y = (0.4 * Rt * sind(P0_angle)) + (0.4 * Rt) + Rt;   %finding y value of P0

%Point P2 (last point on parabola - end of nozzle)
P2_x = eighty_percent_bell_contour_length;              %finding x value of P2
P2_y = exit_radius;                                     %finding y value of P2

%Point P1 (intersection of tangent lines of P0 and P2)
slope_parabola_start = tand(initial_parabola_angle);         %slope at point where parabola portion begins
slope_parabola_end = tand(final_parabola_angle) ;            %slope at point where parabola portion ends
slope_intercept_1 = P0_y - (slope_parabola_start * P0_x);   %intercept of straight line from point P0 with slope at the point
slope_intercept_2 = P2_y - (slope_parabola_end * P2_x);     %intercept of straight line from point P2 with slope at the point
P1_x = (slope_intercept_2 - slope_intercept_1) / (slope_parabola_start - slope_parabola_end); %finding x value of P1
P1_y = ((slope_parabola_start * slope_intercept_2) - (slope_parabola_end * slope_intercept_1)) / (slope_parabola_start - slope_parabola_end); %finding y value of P1

%DIVERGING PARABOLA PORTION
t = linspace(0,1);                                                      %when t = 1, last point on parabola is (length of nozzle, exit radius)
x_parabola = (((1-t).^2).*P0_x) + (2.*(1-t).*t.*P1_x) + ((t.^2).*P2_x); %Quadratic Bezier Curve x portion
y_parabola = (((1-t).^2).*P0_y) + (2.*(1-t).*t.*P1_y) + ((t.^2).*P2_y); %Quadratic Bezier Curve y portion

%CONVERGING PARABOLA SECTION
%Determining parameters for Quadratic Bezier Curve
%Point C0 (point of inflection between combustion converging portion and 1.5 radius)
C0_x = min(x_1_5);      %finding x value of C0
C0_y = max(y_1_5);      %finding y value of C0

%Point C2 (end of combustion chamber (start of convergence))
C2_x = min(x_1_5) - length_chamber_convergence_section;     %finding x value of C2
C2_y = chamber_radius;                                      %finding y value of C2

%Point C1 (intersection of tangent lines of C0 and C2)
Cslope_parabola_start = tand(-135+90);                           %slope at point where parabola portion begins
Cslope_parabola_end = tand(-180);                                %slope at point where parabola portion ends
Cslope_intercept_1 = C0_y - (Cslope_parabola_start * C0_x);     %intercept of straight line from point C0 with slope at the point
Cslope_intercept_2 = C2_y - (Cslope_parabola_end * C2_x);       %intercept of straight line from point C2 with slope at the point

C1_x = (Cslope_intercept_2 - Cslope_intercept_1) / (Cslope_parabola_start - Cslope_parabola_end); %finding x value of C1
C1_y = ((Cslope_parabola_start * Cslope_intercept_2) - (Cslope_parabola_end * Cslope_intercept_1)) / (Cslope_parabola_start - Cslope_parabola_end); %finding y value of C1

%PARABOLA PORTION
t = linspace(0,1); %when t = 1, last point on parabola has y value of combustion chamber radius
Cx_parabola = (((1-t).^2).*C0_x) + (2.*(1-t).*t.*C1_x) + ((t.^2).*C2_x); %Quadratic Bezier Curve x portion
Cy_parabola = (((1-t).^2).*C0_y) + (2.*(1-t).*t.*C1_y) + ((t.^2).*C2_y); %Quadratic Bezier Curve y portion

%CHAMBER PORTION (defined as not including converging section)
x_chamber = linspace((C2_x - Lc), C2_x);              %creating x values of chamber line
y_chamber = linspace(chamber_radius, chamber_radius); %creating y values of chamber line
%plot(x_chamber, y_chamber)

%CURVE
x = [x_chamber(1:end-1) fliplr(Cx_parabola) x_1_5(2:end) x_0_4(2:end) x_parabola(2:end)].';
r = [y_chamber(1:end-1) fliplr(Cy_parabola) y_1_5(2:end) y_0_4(2:end) y_parabola(2:end)].';
% Zcoord = zeros(length(Xcoord), 1);

%

% Code below for exporting into solidworks
%curve = [ Xcoord Ycoord Zcoord ];
%dlmwrite('contour.txt',curve,'delimiter','\t','precision',6)

%figure 
%plot(Xcoord,Ycoord, 'r-'); axis equal
end

