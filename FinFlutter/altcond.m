function [ dens, Talt, pres ] = altcond(Vval,alt,g)
M = 0.0289644;                         % Molar Mass of Air
R = 8.31447;                           % Ideal Gas Constant
Va = Vval;                             %  Velocity 
   if alt<50
     L = -0.0065;
     T0 = 288.15;
     p0 = 101; 
    Talt = T0;                                     %Temperature at altitude
    pres = p0;            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;                                 % Find Reynolds Number
    endif
   if 50<alt<11000
     L = -0.0065;
     T0 = 288.15;
     p0 = 101; 
    Talt = T0+L*(alt-0);                                     %Temperature at altitude
    pres = (p0*1000*(T0/(Talt))^((g*M)/(R*L)));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;                                 % Find Reynolds Number
  endif
  if 20000>alt>11000
    L=0;
    T0=216.65;
    p0=23;
    Talt = T0+L*(alt-11000);                                     %Temperature at altitude
    pres = p0*1000*exp((-g*M*(alt+.01-11000))/(R*T0));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  endif
  if 32000>alt>20000
    L=0.001;
    T0=216.65;
    p0=5.5
    Talt = T0+L*(alt-20000);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  endif
  if 47000>alt>32000
    L=0.0028;
    T0=228.65;
    p0=0.868
    Talt = T0+L*(alt-32000);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  endif
  if 51000>alt>47000
    L=0.0;
    T0=270.65;
    p0=0.110
    Talt = T0+L*(alt-47000);                                     %Temperature at altitude
    pres = p0*1000*exp((-g*M*(alt+.01-47000))/(R*T0));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  endif
  if 71000>alt>51000
    L=-0.0028;
    T0=270.65;
    p0=0.067
    Talt = T0+L*(alt-51000);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  endif
  if alt>71000
    L=-0.002;
    T0=214.65;
    p0=0.004;
    Talt = T0+L*(alt-71000);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  endif
endfunction
