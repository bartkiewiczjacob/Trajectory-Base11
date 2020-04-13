clear all
   for z = 1:1:20000
    [T1(z),a1(z),P1(z),rho1(z)]= altcond2(z,298.15,0);
   [T2(z),a2(z),P2(z),rho2(z)]= atmosisa(z);
    z1(z)=z;
   end
   
   figure("Name","Temperature")
   title("Comparing Temperature")
   xlabel('Altitude (m)')
   ylabel('Temperature (K)')
   hold on
   plot(z1,T1)
   plot(z1,T2, '--')
   legend('altcond2','atmosisa')
   hold off
     
   figure("Name","Speed of Sound")
   title("Comparing Speed of Sound")
   xlabel('Altitude (m)')
   ylabel('Speed of Sound (m/s)')
   hold on
   plot(z1,a1)
   plot(z1,a2, '--')
   legend('altcond2','atmosisa')
   hold off
   
   figure("Name","Pressure")
   title("Comparing Pressure")
   xlabel('Altitude (m)')
   ylabel('Pressure (Pa)')
   hold on
   plot(z1,P1)
   plot(z1,P2, '--')
   legend('altcond2','atmosisa')
   hold off
   
   figure("Name","Density")
   title("Comparing Density")
   xlabel('Altitude (m)')
   ylabel('Density (kg/m^3)')
   hold on
   plot(z1,rho1)
   plot(z1,rho2, '--')
   legend('altcond2','atmosisa')
   hold off