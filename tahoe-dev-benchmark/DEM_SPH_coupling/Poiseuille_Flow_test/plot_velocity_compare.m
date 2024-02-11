% this is used to plot the simulation velocity Vx against Z coordinate,
% with comparing to the analytical result for the Poiseuille Flow as in John
% Morris's paper in 1996

clear
close all

data_22500 = load('velocity_2250');
data_22500 = data_22500(1:51,:);
Z_22500 = data_22500(:,6);
Vx_22500 = data_22500(:,7);

data_45000 = load('velocity_4500');
data_45000 = data_45000(1:51,:);
Z_45000 = data_45000(:,6);
Vx_45000 = data_45000(:,7);

data_112500 = load('velocity_11250');
data_112500 = data_112500(1:51,:);
Z_112500 = data_112500(:,6);
Vx_112500 = data_112500(:,7);

data_225000 = load('velocity_22500');
data_225000 = data_225000(1:51,:);
Z_225000 = data_225000(:,6);
Vx_225000 = data_225000(:,7);

data_1000000 = load('velocity_100000');
data_1000000 = data_1000000(1:51,:);
Z_1000000 = data_1000000(:,6);
Vx_1000000 = data_1000000(:,7);

% analytical solution
nu = 10^(-6);
L = 0.001;
rho = 1000;
F = 10^(-4);

Z_ana = (0:L/200:L);

syms n
for ii=1:length(Z_ana)
   
   ii 
    
   t = 0.0225; 
   Vx_22500_ana(ii) = F/(2*nu)*Z_ana(ii)*(Z_ana(ii)-L) + symsum(4*F*L*L/(nu*pi^3*(2*n+1)^3)*sin(pi*Z_ana(ii)/L*(2*n+1))*exp(-(2*n+1)^2*pi^2*nu/(L*L)*t), 0, 20);
   
   t = 0.045; 
   Vx_45000_ana(ii) = F/(2*nu)*Z_ana(ii)*(Z_ana(ii)-L) + symsum(4*F*L*L/(nu*pi^3*(2*n+1)^3)*sin(pi*Z_ana(ii)/L*(2*n+1))*exp(-(2*n+1)^2*pi^2*nu/(L*L)*t), 0, 20);
   
   t = 0.1125; 
   Vx_112500_ana(ii) = F/(2*nu)*Z_ana(ii)*(Z_ana(ii)-L) + symsum(4*F*L*L/(nu*pi^3*(2*n+1)^3)*sin(pi*Z_ana(ii)/L*(2*n+1))*exp(-(2*n+1)^2*pi^2*nu/(L*L)*t), 0, 20);
   
   t = 0.225; 
   Vx_225000_ana(ii) = F/(2*nu)*Z_ana(ii)*(Z_ana(ii)-L) + symsum(4*F*L*L/(nu*pi^3*(2*n+1)^3)*sin(pi*Z_ana(ii)/L*(2*n+1))*exp(-(2*n+1)^2*pi^2*nu/(L*L)*t), 0, 20);
   
   t = 1; 
   Vx_1000000_ana(ii) = F/(2*nu)*Z_ana(ii)*(Z_ana(ii)-L) + symsum(4*F*L*L/(nu*pi^3*(2*n+1)^3)*sin(pi*Z_ana(ii)/L*(2*n+1))*exp(-(2*n+1)^2*pi^2*nu/(L*L)*t), 0, 20);
    
end


figure(1)

grid on
hold on
plot(Z_22500*1000, Vx_22500*10^5, 'k.')
plot(Z_ana*1000, -Vx_22500_ana*10^5, 'k')

plot(Z_45000*1000, Vx_45000*10^5, 'k.')
plot(Z_ana*1000, -Vx_45000_ana*10^5, 'k')

plot(Z_112500*1000, Vx_112500*10^5, 'k.')
plot(Z_ana*1000, -Vx_112500_ana*10^5, 'k')

plot(Z_225000*1000, Vx_225000*10^5, 'k.')
plot(Z_ana*1000, -Vx_225000_ana*10^5, 'k')

plot(Z_1000000*1000, Vx_1000000*10^5, 'k.')
plot(Z_ana*1000, -Vx_1000000_ana*10^5, 'k')

xlabel('z (10^{-3}m)')
ylabel('v_x (10^{-5}m/s)')
legend('SPH','Series Solution')



