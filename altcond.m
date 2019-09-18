function [ RE, dens, Talt ] = altcond(Vval,zt,g)
M = 28.9644;                         % Molar Mass of Air
R = 89494.596;                           % Ideal Gas Constant
Va = Vval;                                        %  Velocity 
   if zt<50
     L = -0.0065*0.3048;
     T0 = 288.15;
     p0 = 101*0.2953; 
    d0 = 1.2250*0.00194032;
    Talt = T0;                                     %Temperature at altitude
    pres = p0;            %Air Pressure at altitude
    dens = d0*(T0/(T0+L*(zt-0)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;                                 % Find Reynolds Number
    end
   if (zt>50)&&(zt<36089)
     L = -0.0065*0.3048;
     T0 = 288.15;
     p0 = 101*0.2953; 
     d0 = 1.2250*0.00194032;
    Talt = T0+L*(zt-0);                                     %Temperature at altitude
    pres = (p0*1000*(T0/(Talt))^((g*M)/(R*L)));            %Air Pressure at altitude
   dens = d0*(T0/(T0+L*(zt-0)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;                                 % Find Reynolds Number
  end
  if (36089>zt)&&(zt>65617)
    L=0;
    T0=216.65;
    p0=23*0.2953;
    d0=0.36391*0.00194032;
    Talt = T0+L*(zt-36089);                                     %Temperature at altitude
    pres = p0*1000*exp((-g*M*(zt+.01-36089))/(R*T0));            %Air Pressure at altitude
    dens = d0*exp(-g*M*(zt-36089)/(R*T0));                          %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  end
  if (104987>zt)&&(zt>65617)
    L=0.001*0.3048;
    T0=216.65;
    p0=5.5*0.2953;
    d0=0.08803*0.00194032;
    Talt = T0+L*(zt-65617);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = d0*(T0/(T0+L*(zt-65617)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  end
  if (154199>zt)&&(zt>104987)
    L=0.0028*0.3048;
    T0=228.65;
    p0=0.868*0.2953;
    d0 = 0.01322*0.00194032;
    Talt = T0+L*(zt-104987);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = d0*(T0/(T0+L*(zt-104987)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  end
  if (167323>zt)&&(zt>154199)
    L=0.0;
    T0=270.65;
    p0=0.110*0.2953;
    Talt = T0+L*(zt-154199);                                     %Temperature at altitude
    pres = p0*1000*exp((-g*M*(zt+.01-154199))/(R*T0));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  end
  if (232940>zt)&&(zt>167323)
    L=-0.0028*0.3048;
    T0=270.65;
    p0=0.067*0.2953;
    d0 = 0.00086*0.00194032;
    Talt = T0+L*(zt-167323);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = d0*(T0/(T0+L*(zt-167323)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  end
  if (zt>232940)&&(zt<86000)
    L=-0.002*0.3048;
    T0=214.65;
    p0=0.004*0.2953;
    d0 = 0.000064*0.00194032;
    Talt = T0+L*(zt-232940);                                     %Temperature at altitude
    pres = p0*1000*(T0/(Talt))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = d0*(T0/(T0+L*(zt-232940)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;   
  end
  if zt>86000
    L=0;
    T0=186.87;
    p0 = 0.003734*0.2953; 
    d0 = 0.000064*0.00194032;
    Talt = T0+L*(zt-0);                                     %Temperature at altitude
    pres = (p0*1000*(T0/(Talt))^((g*M)/(R*L)));            %Air Pressure at altitude
    dens = d0*(T0/(T0+L*(zt-0)))^(1+(g*M/(R*L)));                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    RE = (Va*dens*1)/mu;         
  end
end