function [ Talt,Va,pres,dens ] = altcond2(zt,Ti,zi)
%
%   altcond2 - Altitude conditions given an initial temperature
%       and an initial altitude for that reading. It is similar
%       to altcond1 which was lacking documentation.
%
%   INPUTS:
%       zt - Altitude in meters above sea level of the desired
%            conditions
%       Ti - Temperature measured at the initial altitude in K-
%            elvin
%       zi - Initial altitude of the initial temperature in me-
%            ters
%
%   OUTPUTS:
%       Talt - Temperature at altitude, Kelvin
%       Va - Velocity of the speed of sount meters/second
%       pres - Air pressure, kPa
%       dens - Air density, kg/m^3
%
%       BY: Joseph McGow-Russell
%       Ver 1 - April 13, 2020
%
M = 0.0289644;                         % Molar Mass of Air kg/mol
R = 8.31447;                           % Ideal Gas Constant kPa/l 
g=9.81;                                % Gravitational Constant m/s^2
   if 0<=zt && zt<11000
     L = -0.0065;                      % Temperature lapse in this layer.
     T0 = Ti-L*(zi);
     %Calculates the initial temperature for this layer
     p0 = 101325; 
    Talt = T0+L*(zt);                                    % Temperature at altitude
    pres = p0*(T0/((T0+L*zt)))^((g*M)/(R*L));            % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                            % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound
   end
  if 11000<=zt && zt<20000
    L=0;                                % Temperature lapse in this layer.
    T0=Ti-(-0.0065*zi)+(-0.0065*11000);
    %Calculates the initial temperature for this layer
    p0=22632.1;
    Talt = T0+L*(zt-11000);                              % Temperature at altitude
    pres = p0*exp((-g*M*(zt-11000))/(R*T0));             % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                            % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound
  end
  if 20000<=zt && zt<32000
    L=0.001;                           % Temperature lapse in this layer.
    T0=Ti-(-0.0065*zi)+(-0.0065*11000)+(0*20000);
    %Calculates the initial temperature for this layer
    p0=5474.89;
    Talt = T0+L*(zt-20000);                              % Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-20000))))^((g*M)/(R*L));    % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                            % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound 
  end
  if 32000<=zt && zt<47000
    L=0.0028;                          % Temperature lapse in this layer.
    T0=Ti-(-0.0065*zi)+(-0.0065*11000)+(0*(20000-11000))+(0.001*(32000-20000));
    %Calculates the initial temperature for this layer
    p0=868.02;
    Talt = T0+L*(zt-32000);                              % Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-32000))))^((g*M)/(R*L));    % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                            % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound
  end
  if 47000<=zt && zt<51000
    L=0;                                % Temperature lapse in this layer.
    T0=Ti-(-0.0065*zi)+(-0.0065*11000)+(0*(20000-11000))+(0.001*(32000-20000))+(0.0028*(47000-32000));
    %Calculates the initial temperature for this layer
    p0=110.91;
    Talt = T0+L*(zt-47000);                             % Temperature at altitude
    pres = p0*exp((-g*M*(zt-47000))/(R*T0));            % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
   Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound 
  end
  if 51000<=zt && zt<71000
    L=-0.0028;                             % Temperature lapse in this layer
    T0=Ti-(-0.0065*zi)+(-0.0065*11000)+(0*(20000-11000))+(0.001*(32000-20000))+(0.0028*(47000-32000))+(0*(51000-47000));
    %Calculates the initial temperature for this layer
    p0=66.94;
    Talt = T0+L*(zt-51000);                              % Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-51000))))^((g*M)/(R*L));    % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                            % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound
  end
  if 71000<=zt && zt<100000
    L=-0.002;                                % Temperature lapse in this layer
    T0=Ti-(-0.0065*zi)+(-0.0065*11000)+(0*(20000-11000))+(0.001*(32000-20000))+(0.0028*(47000-32000))+(0*(51000-47000))+(-0.0028*(71000-51000));
    %Calculates the initial temperature for this layer
    p0=3.96;
    Talt = T0+L*(zt-71000);                             % Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-71000))))^((g*M)/(R*L));   % Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           % Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
   Va = sqrt(1.4*287.058*Talt);                         % Velocity of Speed of Sound
  end
  if 100000<zt
      %%Assume we are in SPACE%%
      Talt=0;
      pres=0;
      dens=0;
      Va=0;
  end
end
