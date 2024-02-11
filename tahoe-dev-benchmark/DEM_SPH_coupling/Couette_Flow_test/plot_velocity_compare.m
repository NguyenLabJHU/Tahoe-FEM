% this is used to plot the simulation velocity Vx against Z coordinate,
% with comparing to the analytical result for the Couette Flow as in John
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
V0 = 1.25*10^(-5);

Z_ana = (0:L/200:L);

syms n
for ii=1:length(Z_ana)
   
   ii 
    
   t = 0.0225; 
   Vx_22500_ana(ii) = V0/L*Z_ana(ii) + symsum(2*V0/(n*pi)*(-1)^n*sin(n*pi/L*Z_ana(ii))*exp(-nu*n*n*pi*pi/L/L*t), 1, 20);
   
   t = 0.045; 
   Vx_45000_ana(ii) = V0/L*Z_ana(ii) + symsum(2*V0/(n*pi)*(-1)^n*sin(n*pi/L*Z_ana(ii))*exp(-nu*n*n*pi*pi/L/L*t), 1, 20);
   
   t = 0.1125; 
   Vx_112500_ana(ii) = V0/L*Z_ana(ii) + symsum(2*V0/(n*pi)*(-1)^n*sin(n*pi/L*Z_ana(ii))*exp(-nu*n*n*pi*pi/L/L*t), 1, 20);
   
   t = 0.225; 
   Vx_225000_ana(ii) = V0/L*Z_ana(ii) + symsum(2*V0/(n*pi)*(-1)^n*sin(n*pi/L*Z_ana(ii))*exp(-nu*n*n*pi*pi/L/L*t), 1, 20);
   
   t = 1; 
   Vx_1000000_ana(ii) = V0/L*Z_ana(ii) + symsum(2*V0/(n*pi)*(-1)^n*sin(n*pi/L*Z_ana(ii))*exp(-nu*n*n*pi*pi/L/L*t), 1, 20);
    
end


figure(1)

grid on
hold on
plot(Z_22500*1000, Vx_22500*10^5, 'k.')
plot(Z_ana*1000, Vx_22500_ana*10^5, 'k')

plot(Z_45000*1000, Vx_45000*10^5, 'k.')
plot(Z_ana*1000, Vx_45000_ana*10^5, 'k')

plot(Z_112500*1000, Vx_112500*10^5, 'k.')
plot(Z_ana*1000, Vx_112500_ana*10^5, 'k')

plot(Z_225000*1000, Vx_225000*10^5, 'k.')
plot(Z_ana*1000, Vx_225000_ana*10^5, 'k')

plot(Z_1000000*1000, Vx_1000000*10^5, 'k.')
plot(Z_ana*1000, Vx_1000000_ana*10^5, 'k')

xlabel('z (10^{-3}m)')
ylabel('v_x (10^{-5}m/s)')
legend('SPH','Series Solution')



