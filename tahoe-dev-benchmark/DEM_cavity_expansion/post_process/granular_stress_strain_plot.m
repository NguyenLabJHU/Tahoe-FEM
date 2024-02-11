close all
clear all

sim = load('exp_progress_RVE');
timestep = 5.0e-7;
print_inter = 1000;
sim_inter = 3;
sim = sim(1:sim_inter:end,:);
sim = sim(2:1:end,:);

Bagi_linear_11=sim(:,36); %sample epsilon_w
Bagi_linear_12=sim(:,37);
Bagi_linear_13=sim(:,38);
Bagi_linear_21=sim(:,39);
Bagi_linear_22=sim(:,40); %sample epsilon_l
Bagi_linear_23=sim(:,41);
Bagi_linear_31=sim(:,42);
Bagi_linear_32=sim(:,43);
Bagi_linear_33=sim(:,44); %sample epsilon_h
Bagi_linear_v=(Bagi_linear_11+Bagi_linear_22+Bagi_linear_33);

old_epsilon_w=-sim(:,32); %sample epsilon_w
old_epsilon_l=-sim(:,33); %sample epsilon_l
old_epsilon_h=-sim(:,34); %sample epsilon_h
old_epsilon_v=-sim(:,35);


% % % finite granular strain
finite_11 = sim(:,51);
finite_12 = sim(:,52);
finite_13 = sim(:,53);
finite_21 = sim(:,54);
finite_22 = sim(:,55);
finite_23 = sim(:,56);
finite_31 = sim(:,57);
finite_32 = sim(:,58);
finite_33 = sim(:,59);
finite_v = finite_11+finite_22+finite_33;

% mixed finite granular strain
lag_11 = sim(:,60);
lag_12 = sim(:,61);
lag_13 = sim(:,62);
lag_21 = sim(:,63);
lag_22 = sim(:,64);
lag_23 = sim(:,65);
lag_31 = sim(:,66);
lag_32 = sim(:,67);
lag_33 = sim(:,68);
lag_v = lag_11+lag_22+lag_33;

% eulerian granular strain
eulerian_11 = sim(:,69);
eulerian_12 = sim(:,70);
eulerian_13 = sim(:,71);
eulerian_21 = sim(:,72);
eulerian_22 = sim(:,73);
eulerian_23 = sim(:,74);
eulerian_31 = sim(:,75);
eulerian_32 = sim(:,76);
eulerian_33 = sim(:,77);
eulerian_v = eulerian_11+eulerian_22+eulerian_33;

% mixed finite granular strain
euler_11 = sim(:,78);
euler_12 = sim(:,79);
euler_13 = sim(:,80);
euler_21 = sim(:,81);
euler_22 = sim(:,82);
euler_23 = sim(:,83);
euler_31 = sim(:,84);
euler_32 = sim(:,85);
euler_33 = sim(:,86);
euler_v = euler_11+euler_22+euler_33;

% % logrithmic strain
% log_epsilon_w=sim(:,87); %sample epsilon_w
% log_epsilon_l=sim(:,88); %sample epsilon_l
% log_epsilon_h=sim(:,89); %sample epsilon_h
% log_epsilon_v=sim(:,90);

% Bagi dudx
Bagi_dudx_11 = sim(:,91);
Bagi_dudx_12 = sim(:,92);
Bagi_dudx_13 = sim(:,93);
Bagi_dudx_21 = sim(:,94);
Bagi_dudx_22 = sim(:,95);
Bagi_dudx_23 = sim(:,96);
Bagi_dudx_31 = sim(:,97);
Bagi_dudx_32 = sim(:,98);
Bagi_dudx_33 = sim(:,99);

% Lagrangian dudx
Lag_dudx_11 = sim(:,100);
Lag_dudx_12 = sim(:,101);
Lag_dudx_13 = sim(:,102);
Lag_dudx_21 = sim(:,103);
Lag_dudx_22 = sim(:,104);
Lag_dudx_23 = sim(:,105);
Lag_dudx_31 = sim(:,106);
Lag_dudx_32 = sim(:,107);
Lag_dudx_33 = sim(:,108);

% Eulerian dudx
Euler_dudx_11 = sim(:,109);
Euler_dudx_12 = sim(:,110);
Euler_dudx_13 = sim(:,111);
Euler_dudx_21 = sim(:,112);
Euler_dudx_22 = sim(:,113);
Euler_dudx_23 = sim(:,114);
Euler_dudx_31 = sim(:,115);
Euler_dudx_32 = sim(:,116);
Euler_dudx_33 = sim(:,117);

% strain based on spatial deformation rate tensor, July 9, 2013
rate_11 = sim(:,127)*1000;
rate_12 = sim(:,128)*1000;
rate_13 = sim(:,129)*1000;
rate_21 = sim(:,130)*1000;
rate_22 = sim(:,131)*1000;
rate_23 = sim(:,132)*1000;
rate_31 = sim(:,133)*1000;
rate_32 = sim(:,134)*1000;
rate_33 = sim(:,135)*1000;

rate_strain_v = rate_11+rate_22+rate_33;

% % above is for triaxial compression results


% calculate Bagi's finite strain based Bagi_dudx
iden = [1 0 0; 0 1 0; 0 0 1];   % for mapped Lagrangian
n_steps = length(Bagi_dudx_11);
time_vec = [1:n_steps]*print_inter*sim_inter*timestep;
for i=1:n_steps
Bagi_dudx = [Bagi_dudx_11(i) Bagi_dudx_12(i) Bagi_dudx_13(i);
             Bagi_dudx_21(i) Bagi_dudx_22(i) Bagi_dudx_23(i);
             Bagi_dudx_31(i) Bagi_dudx_32(i) Bagi_dudx_33(i)];
Bagi_finite = (Bagi_dudx+Bagi_dudx'-Bagi_dudx'*Bagi_dudx)/2;
Bagi_11(i) = Bagi_finite(1,1);
Bagi_22(i) = Bagi_finite(2,2);
Bagi_33(i) = Bagi_finite(3,3);

Bagi_12(i) = Bagi_finite(1,2);
Bagi_13(i) = Bagi_finite(1,3);
Bagi_23(i) = Bagi_finite(2,3);

Bagi_v(i) = Bagi_finite(1,1)+Bagi_finite(2,2)+Bagi_finite(3,3);


% check to see why Bagi_linear == Euler
Bagi_linear = (Bagi_dudx+Bagi_dudx')/2;
Bagi_linear_11(i) = Bagi_linear(1,1);
Bagi_linear_22(i) = Bagi_linear(2,2);
Bagi_linear_33(i) = Bagi_linear(3,3);

Bagi_linear_12(i) = Bagi_linear(1,2);
Bagi_linear_13(i) = Bagi_linear(1,3);
Bagi_linear_23(i) = Bagi_linear(2,3);

Bagi_linear_v(i) = Bagi_linear(1,1)+Bagi_linear(2,2)+Bagi_linear(3,3);





% finite Lagrangian strain based Lag_dudx
Lag_dudx = [Lag_dudx_11(i) Lag_dudx_12(i) Lag_dudx_13(i);
            Lag_dudx_21(i) Lag_dudx_22(i) Lag_dudx_23(i);
            Lag_dudx_31(i) Lag_dudx_32(i) Lag_dudx_33(i)];
Lag_finite = (Lag_dudx+Lag_dudx'+Lag_dudx'*Lag_dudx)/2;
% finite_11_dudx(i) = Lag_finite(1,1);
% finite_22_dudx(i) = Lag_finite(2,2);
% finite_33_dudx(i) = Lag_finite(3,3);
% 
% finite_12_dudx(i) = Lag_finite(1,2);


Lag_11(i) = Lag_finite(1,1);
Lag_22(i) = Lag_finite(2,2);
Lag_33(i) = Lag_finite(3,3);

Lag_12(i) = Lag_finite(1,2);
Lag_13(i) = Lag_finite(1,3);
Lag_23(i) = Lag_finite(2,3);

Lag_v(i) = Lag_finite(1,1)+Lag_finite(2,2)+Lag_finite(3,3);

Lag_linear = (Lag_dudx+Lag_dudx')/2;
Lag_linear_11(i) = Lag_linear(1,1);
Lag_linear_22(i) = Lag_linear(2,2);
Lag_linear_33(i) = Lag_linear(3,3);

Lag_linear_12(i) = Lag_linear(1,2);
Lag_linear_13(i) = Lag_linear(1,3);
Lag_linear_23(i) = Lag_linear(2,3);

Lag_linear_v(i) = Lag_linear(1,1)+Lag_linear(2,2)+Lag_linear(3,3);



Lag_F = Lag_dudx+iden;  % for mapped Lagrangian

% Lag_HOT = Lag_dudx'*Lag_dudx/2;
% Lag_HOT_11(i) = Lag_HOT(1,1);
% Lag_HOT_22(i) = Lag_HOT(2,2);
% Lag_HOT_33(i) = Lag_HOT(3,3);
% 
% Lag_HOT_12(i) = Lag_HOT(1,2);
% Lag_HOT_13(i) = Lag_HOT(1,3);
% Lag_HOT_23(i) = Lag_HOT(2,3);
% 
% Lag_HOT_v(i) = Lag_HOT(1,1)+Lag_HOT(2,2)+Lag_HOT(3,3);
% finite_v_dudx(i) = Lag_finite(1,1)+Lag_finite(2,2)+Lag_finite(3,3);

% finite Eulerian strain based Euler_dudx
Euler_dudx = [Euler_dudx_11(i) Euler_dudx_12(i) Euler_dudx_13(i);
              Euler_dudx_21(i) Euler_dudx_22(i) Euler_dudx_23(i);
              Euler_dudx_31(i) Euler_dudx_32(i) Euler_dudx_33(i)];
Euler_F = (iden-Euler_dudx)^(-1);
Euler_finite = (Euler_dudx+Euler_dudx'-Euler_dudx'*Euler_dudx)/2;% for mapped Lagrangian
map_Lag = Lag_F'*Euler_finite*Lag_F;
map_Lag_11(i) = map_Lag(1,1);
map_Lag_22(i) = map_Lag(2,2);
map_Lag_33(i) = map_Lag(3,3);

map_Lag_12(i) = map_Lag(1,2);
map_Lag_13(i) = map_Lag(1,3);
map_Lag_23(i) = map_Lag(2,3);

map_Lag_v(i) = map_Lag(1,1)+map_Lag(2,2)+map_Lag(3,3);


Euler_11(i) = Euler_finite(1,1);
Euler_22(i) = Euler_finite(2,2);
Euler_33(i) = Euler_finite(3,3);

Euler_12(i) = Euler_finite(1,2);
Euler_13(i) = Euler_finite(1,3);
Euler_23(i) = Euler_finite(2,3);

Euler_v(i) = Euler_finite(1,1)+Euler_finite(2,2)+Euler_finite(3,3);

Euler_linear = (Euler_dudx+Euler_dudx')/2;
Euler_linear_11(i) = Euler_linear(1,1);
Euler_linear_22(i) = Euler_linear(2,2);
Euler_linear_33(i) = Euler_linear(3,3);

Euler_linear_12(i) = Euler_linear(1,2);
Euler_linear_13(i) = Euler_linear(1,3);
Euler_linear_23(i) = Euler_linear(2,3);

Euler_linear_v(i) = Euler_linear(1,1)+Euler_linear(2,2)+Euler_linear(3,3);




% Euler_HOT = -Euler_dudx'*Euler_dudx/2;
% Euler_HOT_11(i) = Euler_HOT(1,1);
% Euler_HOT_22(i) = Euler_HOT(2,2);
% Euler_HOT_33(i) = Euler_HOT(3,3);
% 
% Euler_HOT_12(i) = Euler_HOT(1,2);
% Euler_HOT_13(i) = Euler_HOT(1,3);
% Euler_HOT_23(i) = Euler_HOT(2,3);
% 
% Euler_HOT_v(i) = Euler_HOT(1,1)+Euler_HOT(2,2)+Euler_HOT(3,3);
% eulerian_v_dudx(i) = Euler_finite(1,1)+Euler_finite(2,2)+Euler_finite(3,3);



end

% Euler_linear_11 = eulerian_11-Euler_HOT_11';
% Euler_linear_22 = eulerian_11-Euler_HOT_22';
% Euler_linear_33 = eulerian_11-Euler_HOT_33';
% Euler_linear_12 = eulerian_11-Euler_HOT_12';
% Euler_linear_13 = eulerian_11-Euler_HOT_13';
% Euler_linear_23 = eulerian_11-Euler_HOT_23';

% Lag_linear_11 = finite_11-Lag_HOT_11';
% Lag_linear_22 = finite_11-Lag_HOT_22';
% Lag_linear_33 = finite_11-Lag_HOT_33';
% Lag_linear_12 = finite_11-Lag_HOT_12';
% Lag_linear_13 = finite_11-Lag_HOT_13';
% Lag_linear_23 = finite_11-Lag_HOT_23';

% calculate octahedral shear strain
% Lagrangian finite 
Lag_oct = sqrt(1/9*((Lag_11-Lag_22).^2+(Lag_22-Lag_33).^2+(Lag_11-Lag_33).^2)+...
          2/3*(Lag_12.^2+Lag_13.^2+Lag_23.^2));
      
% linear Lagrangian 
Lag_linear_oct = sqrt(1/9*((Lag_linear_11-Lag_linear_22).^2+(Lag_linear_22-Lag_linear_33).^2+(Lag_linear_11-Lag_linear_33).^2)+...
                 2/3*(Lag_linear_12.^2+Lag_linear_13.^2+Lag_linear_23.^2));

% Eulerian finite 
Euler_oct = sqrt(1/9*((Euler_11-Euler_22).^2+(Euler_22-Euler_33).^2+(Euler_11-Euler_33).^2)+...
            2/3*(Euler_12.^2+Euler_13.^2+Euler_23.^2));
      
% linear Eulerian
Euler_linear_oct = sqrt(1/9*((Euler_linear_11-Euler_linear_22).^2+(Euler_linear_22-Euler_linear_33).^2+(Euler_linear_11-Euler_linear_33).^2)+...
                   2/3*(Euler_linear_12.^2+Euler_linear_13.^2+Euler_linear_23.^2));
             
% Bagi finite 
Bagi_oct = sqrt(1/9*((Bagi_11-Bagi_22).^2+(Bagi_22-Bagi_33).^2+(Bagi_11-Bagi_33).^2)+...
           2/3*(Bagi_12.^2+Bagi_13.^2+Bagi_23.^2));
      
% linear Bagi
Bagi_linear_oct = sqrt(1/9*((Bagi_linear_11-Bagi_linear_22).^2+(Bagi_linear_22-Bagi_linear_33).^2+(Bagi_linear_11-Bagi_linear_33).^2)+...
                 2/3*(Bagi_linear_12.^2+Bagi_linear_13.^2+Bagi_linear_23.^2));
             
% map Lagrangian finite 
map_Lag_oct = sqrt(1/9*((map_Lag_11-map_Lag_22).^2+(map_Lag_22-map_Lag_33).^2+(map_Lag_11-map_Lag_33).^2)+...
              2/3*(map_Lag_12.^2+map_Lag_13.^2+map_Lag_23.^2));
          
% strain based on spatial deformation rate tensor
rate_strain_oct = sqrt(1/9*((rate_11-rate_22).^2+(rate_22-rate_33).^2+(rate_11-rate_33).^2)+...
          2/3*(rate_12.^2+rate_13.^2+rate_23.^2));  
          
               
% 
% for i=1:timesteps
%     
%    e_Bagi_epsilon_11(i) = new_epsilon_11/log_epsilon_l-1;
%    e_li_Lag_11(i) = (finite_11-lag_11)/log_epsilon_l(i)-1;
%    e_li_Eul_11(i) = (eulerian_11-euler_11)/log_epsilon_l(i)-1;
%     
%     
% end

% plot volumetric strain
figure(1)

hold on
grid on
plot(time_vec, Bagi_v,'ro');
plot(time_vec, Bagi_linear_v,'rx');

plot(time_vec, Lag_v,'bo');
plot(time_vec, Lag_linear_v,'bx');

plot(time_vec, Euler_v,'go');
% plot(time_vec, eulerian_v-Euler_HOT_v','gx');
plot(time_vec, Euler_linear_v,'gx');

% plot(time_vec, map_Lag_v, 'd--')
plot(time_vec, rate_strain_v, 'rd')
xlabel('time (s)')
ylabel('granular strain')
title('volumetric strain')
legend('Bagi finite e_v', 'linear Bagi e_v',...
       'Lagrangian E_v', 'linear Lagrangian E_v', ...
       'Eulerian e_v', 'linear Eulerian e_v', ...
       'rate strain e_v');

   
% plot octahedral shear strain
figure(2)

hold on
grid on
plot(time_vec, Bagi_oct,'ro');
plot(time_vec, Bagi_linear_oct,'rx');

plot(time_vec, Lag_oct,'bo');
plot(time_vec, Lag_linear_oct,'bx');

plot(time_vec, Euler_oct,'go');
plot(time_vec, Euler_linear_oct,'gx');

% plot(time_vec, map_Lag_oct, 'd--')
plot(time_vec, rate_strain_oct, 'rd')
xlabel('time (s)')
ylabel('granular strain')
title('octahedral shear strain')
legend('Bagi finite oct', 'linear Bagi oct',...
       'Lagrangian oct', 'linear Lagrangian oct', ...
       'Eulerian oct', 'linear Eulerian oct', ...
       'rate strain oct');
   
   
   
   
   
% stress read, calculation and plot
granular_sigma_11=sim(:,18); %sample sigma1_1
granular_sigma_12=sim(:,19); %sample sigma1_2
granular_sigma_13=sim(:,20); %sample sigma2_1
granular_sigma_21=sim(:,21); %sample sigma2_2
granular_sigma_22=sim(:,22); %sample sigma3_1
granular_sigma_23=sim(:,23); %sample sigma3_2
granular_sigma_31=sim(:,24);
granular_sigma_32=sim(:,25);
granular_sigma_33=sim(:,26);
   
mean_stress = (granular_sigma_11+granular_sigma_22+granular_sigma_33)/3;

% octahedral shear stress
stress_oct = sqrt(1/9*((granular_sigma_11-granular_sigma_22).^2+(granular_sigma_22-granular_sigma_33).^2+(granular_sigma_11-granular_sigma_33).^2)+...
             2/3*(granular_sigma_12.^2+granular_sigma_13.^2+granular_sigma_23.^2));

% plot
figure(21)
hold on 
grid on
plot(time_vec(1:end), mean_stress(1:end), 'd--','LineWidth',2);
xlabel('time (s)')
ylabel('granular stress (Pa)')
title('mean stress')
legend('mean stress')

figure(22)
hold on
grid on
plot(time_vec(1:end), stress_oct(1:end), 'd--','LineWidth',2);
xlabel('time (s)')
ylabel('grnaular stress')
title('octahedral stress')
legend('octahedral stress')




% plot stress-strain relation
% plot mean stress - volumetric strain
figure(31)

hold on
grid on
% plot(Bagi_v(1:end), mean_stress(1:end), 'ro--');
% plot(Bagi_linear_v(1:end), mean_stress(1:end), 'rx--');
% 
% plot(Lag_v(1:end), mean_stress(1:end), 'bo--');
% plot(Lag_linear_v(1:end), mean_stress(1:end), 'bx--');

plot(Euler_v(1:end), mean_stress(1:end), 'go--','LineWidth',2);
% plot(Euler_linear_v(1:end), mean_stress(1:end), 'gx--');
plot(rate_strain_v(1:end), mean_stress(1:end), 'rd--', 'LineWidth',2)
xlabel('granular strain')
ylabel('granular stress (Pa)')
title('mean stress vs. volumetric strain')
legend('Eulerian granular strain', 'granular strain by rate form')


   
% plot octahedral shear stress - octahedral shear strain
figure(32)

hold on
grid on
% plot(Bagi_oct(1:end), stress_oct(1:end), 'ro--');
% plot(Bagi_linear_oct(1:end), stress_oct(1:end), 'rx--');
% 
% plot(Lag_oct(1:end), stress_oct(1:end), 'bo--');
% plot(Lag_linear_oct(1:end), stress_oct(1:end), 'bx--');

plot(Euler_oct(1:end), stress_oct(1:end), 'go--','LineWidth',2);
% plot(Euler_linear_oct(1:end), stress_oct(1:end), 'gx--');
plot(rate_strain_oct(1:end), stress_oct(1:end), 'rd--','LineWidth',2);
xlabel('granular strain')
ylabel('granualr stress (Pa)')
title('octahedral shear stress vs. octahedral shear strain')
legend('Eulerian granular strain', 'granular strain by rate form')
% legend('Bagi finite \gamma', 'linear Bagi oct',...
%        'Lagrangian oct', 'linear Lagrangian oct', ...
%        'Eulerian oct', 'linear Eulerian oct');
   
  
   
% spatial dvdx
spatial_dvdx_11 = sim(:,118);
spatial_dvdx_12 = sim(:,119);
spatial_dvdx_13 = sim(:,120);
spatial_dvdx_21 = sim(:,121);
spatial_dvdx_22 = sim(:,122);
spatial_dvdx_23 = sim(:,123);
spatial_dvdx_31 = sim(:,124);
spatial_dvdx_32 = sim(:,125);
spatial_dvdx_33 = sim(:,126);   

for i=1:n_steps
spatial_dvdx = [spatial_dvdx_11(i) spatial_dvdx_12(i) spatial_dvdx_13(i);
                spatial_dvdx_21(i) spatial_dvdx_22(i) spatial_dvdx_23(i);
                spatial_dvdx_31(i) spatial_dvdx_32(i) spatial_dvdx_33(i)];
spatial_rate = (spatial_dvdx+spatial_dvdx')/2;
spatial_rate_11(i) = spatial_rate(1,1);
spatial_rate_22(i) = spatial_rate(2,2);
spatial_rate_33(i) = spatial_rate(3,3);

spatial_rate_12(i) = spatial_rate(1,2);
spatial_rate_13(i) = spatial_rate(1,3);
spatial_rate_23(i) = spatial_rate(2,3);

spatial_rate_v(i) = spatial_rate(1,1)+spatial_rate(2,2)+spatial_rate(3,3);

end
   

rate_oct = sqrt(1/9*((spatial_rate_11-spatial_rate_22).^2+(spatial_rate_22-spatial_rate_33).^2+(spatial_rate_11-spatial_rate_33).^2)+...
           2/3*(spatial_rate_12.^2+spatial_rate_13.^2+spatial_rate_23.^2));
   
% plot
figure(41)
hold on 
grid on
plot(time_vec, spatial_rate_v, 'd--','LineWidth',2);
xlabel('time (s)')
ylabel('deformation rate (s^{-1})')
title('volumetric granular deformation rate')
% legend('volumetric deformation rate')

figure(42)
hold on
grid on
plot(time_vec, rate_oct, 'd--','LineWidth',2);
xlabel('time (s)')
ylabel('deformation rate (s^{-1})')
title('octahedral shear deformation rate')
% legend('octahedral strain rate')   

figure(43)
hold on
grid on
plot(time_vec, spatial_rate_33, 'd--','LineWidth',2)
xlabel('time (s)')
ylabel('strain rate')
title('normal strain rate in z direction')
legend('strain rate 33')
%    

