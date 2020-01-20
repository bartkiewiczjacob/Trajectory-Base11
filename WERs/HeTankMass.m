function [ massTank ] = HeTankMass( pressure , temp , heMass )
% pressure - pressure inside helium tank, PSIA
% temp - temperature of helium in tank, *R
% heMass - mass of helium to be contained in tank, lbm
%
% estimates tank mass built of Al 6061 with hemisherical bulkheads
% 
%% Constants
R = 386.047 ; % Specific gas constant, (ft*lbf)/(lbm*R)
SF = 2 ; % Safety Factor
sigmaYield = 40000 ; % Yield stress , Al 6061-T6
radiusOuter = 3 ; % Maximum outer diameter, in
rho = 0.0975 ; % Density of material, lbm/in^3

%% Math
sigmaMax = sigmaYield / SF ; % Maximum acceptable hoop stress
volume = heMass * R * temp / ( 144 * pressure ) ; % volume of tank, ft^3
volume = volume * ( 12 ^ 3 ) ; % volume of tank, in^3
radii = ( radiusOuter - 0.05 ) : -0.05 : 0 ;
thickness = 0 ;
radiusInner = 0 ;
for index = 1 : length( radii )
    r = radii( index ) ;
    t = radiusOuter - r ;
    radiusEffective = sqrt( ( radiusOuter ^ 2 + r ^ 2 ) / 2 ) ;
    stress = pressure * radiusEffective / t ;
    if stress < sigmaMax 
        thickness = t ;
        radiusInner = r ;
        break ;
    end
end
volumeCylinder = volume - ( 4 / 3 * pi * ( radiusInner ^ 3 ) ) ;
heightCylinder = volumeCylinder / ( pi * ( radiusInner ^ 2 ) ) ;
thicknessHemispheres = thickness / 2 ;
massHemispheres = rho * ( 4 / 3 ) * pi * ( ( radiusInner + thicknessHemispheres ) ^ 3 - radiusInner ^ 3 ) ;
massCylinder = rho * heightCylinder * ( pi * ( radiusOuter ^ 2 - radiusInner ^ 2 ) ) ;
massTank = massHemispheres + massCylinder ;
fprintf( "\nInner Radius = %f in" , radiusInner ) ;
fprintf( "\nOuter Radius = %f in" , radiusOuter ) ;
fprintf( "\nWall Thickness = %f in " , thickness ) ;
fprintf( "\nCylinder Height = %f in" , heightCylinder ) ;
fprintf( "\nTank Mass = %f lbm" , massTank ) ;

        
    
    
