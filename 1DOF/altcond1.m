function [ Talt,Va,dens ] = altcond1(zt)
M = 0.0289644;                         % Molar Mass of Air
R = 8.31447;                           % Ideal Gas Constant
   if zt<20
     L = -0.0065;
     T0 = 288.15;
     p0 = 101; 
    Talt = T0;                                     %Temperature at altitude
    pres = p0;            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                                        %  Velocity 
    RE = (Va*dens*1)/mu;                                 % Find Reynolds Number
    endif
   if 0<zt<11000
     L = -0.0065;
     T0 = 288.15;
     p0 = 101; 
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*abs(zt)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                                        %  Velocity 
  endif
  if 20000>zt>11000
    L=0;
    T0=216.65;
    p0=23;
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*(zt+.1-11000)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                                        %  Velocity 
  endif
  if 32000>zt>20000
    L=0.001;
    T0=216.65;
    p0=5.5;
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*(zt+.1-20000)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt); 
  endif
  if 47000>zt>32000
    L=0.0028;
    T0=228.65;
    p0=0.868;
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*(zt+.1-0)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);
  endif
  if 51000>zt>47000
    L=0.0;
    T0=270.65;
    p0=0.110;
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*(zt+.1-0)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
   Va = sqrt(1.4*287.058*Talt); 
  endif
  if 71000>zt>51000
    L=-0.0028;
    T0=270.65;
    p0=0.067;
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*(zt+.1-0)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);
  endif
  if zt>71000
    L=-0.002;
    T0=214.65;
    p0=0.004;
    Talt = T0+L*zt;                                     %Temperature at altitude
    pres = p0*1000*(T0/((T0+L)*(zt+.1-0)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
   Va = sqrt(1.4*287.058*Talt);
  endif
endfunction
