clc;
wetted_area = pi*1.5*26;
cross_section = pi*(1.5/2)^2;
l_over_d = 26/1.5;
cbd(:, 1) = basedrag(reynolds(:,1), wetted_area, cross_section, l_over_d);
figure(1)
plot(altitude, cbd(:, 1));
title('Base Drag Coeff vs. Altitude');
xlabel('Altitude [ft]');
ylabel('Base Drag Coefficient');
hold on
wetted_area = pi*1*26;
cross_section = pi*(1/2)^2;
l_over_d = 26/1;
cbd(:, 2) = basedrag(reynolds(:,1), wetted_area, cross_section, l_over_d);
plot(altitude, cbd(:, 2));
wetted_area = pi*2*26;
cross_section = pi*(2/2)^2;
l_over_d = 26/2;
cbd(:, 3) = basedrag(reynolds(:,1), wetted_area, cross_section, l_over_d);
plot(altitude, cbd(:, 3));
legend('18 in', '12 in', '24 in');
hold off
grid on;
figure(2)
hold on
plot(altitude, cbd(:, 1).*q);
plot(altitude, cbd(:, 2).*q);
plot(altitude, cbd(:, 3).*q);
title('Base Drag Force vs. Altitude')
xlabel('Altitude [ft]');
ylabel('Base Drag Force [lbf]');
legend('18 in', '12 in', '24 in');
hold off;
grid on



