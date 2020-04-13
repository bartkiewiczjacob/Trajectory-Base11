function [ Talt,Va,pres,dens ] = altcond1(zt)
M = 0.0289644;                         % Molar Mass of Air
R = 8.31447;                           % Ideal Gas Constant
g=9.81;
   if 0<=zt && zt<11000
     L = -0.0065;
     T0 = 288.15;
     p0 = 101325; 
    Talt = T0+L*(zt);                                     %Temperature at altitude
    pres = p0*(T0/((T0+L*zt)))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                                        %  Velocity 
   end
  if 11000<=zt && zt<20000
    L=0;
    T0=216.65;
    p0=22632.1;
    Talt = T0+L*(zt-11000);                                     %Temperature at altitude
    pres = p0*exp((-g*M*(zt-11000))/(R*T0));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);                                        %  Velocity 
  end
  if 20000<=zt && zt<32000
    L=0.001;
    T0=216.65;
    p0=5474.89;
    Talt = T0+L*(zt-20000);                                     %Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-20000))))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt); 
  end
  if 32000<=zt && zt<47000
    L=0.0028;
    T0=228.65;
    p0=868.02;
    Talt = T0+L*(zt-32000);                                     %Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-32000))))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);
  end
  if 47000<=zt && zt<51000
    L=0.0;
    T0=270.65;
    p0=110.91;
    Talt = T0+L*(zt-47000);                                     %Temperature at altitude
    pres = p0*exp((-g*M*(zt-47000))/(R*T0));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
   Va = sqrt(1.4*287.058*Talt); 
  end
  if 51000<=zt && zt<71000
    L=-0.0028;
    T0=270.65;
    p0=66.94;
    Talt = T0+L*(zt-51000);                                    %Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-51000))))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
    Va = sqrt(1.4*287.058*Talt);
  end
  if 71000<=zt && zt<100000
    L=-0.002;
    T0=214.65;
    p0=3.96;
    Talt = T0+L*(zt-71000);                                     %Temperature at altitude
    pres = p0*(T0/((T0+L*(zt-71000))))^((g*M)/(R*L));            %Air Pressure at altitude
    dens = (pres*M)/(R*Talt);                           %Air Density at altitude
    mu = 1.512041288*(Talt^(3/2))/(Talt+120);
   Va = sqrt(1.4*287.058*Talt);
  end
  if 100000<zt
      Talt=0;
      pres=0;
      dens=0;
      Va=0;
  end
end
