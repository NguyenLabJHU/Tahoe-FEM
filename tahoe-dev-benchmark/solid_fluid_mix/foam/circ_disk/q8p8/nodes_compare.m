
%close all
clear all

format short e

%scale factors
scale_stress=1e0; %keep as MPa
scale_stress=1e3; %convert MPa to kPa

%%%%%%%%%%%%%%%%%%%%%%%%%read nodal abaqus data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('disk_axisymm_abaqus_elastic_S22.txt');
time_drained_air_abaqus=(1e3)*sim(:,1); %convert to ms
S22_drained_air_abaqus=scale_stress*(1e-6)*sim(:,2); %convert to MPa
sim=load('disk_axisymm_abaqus_poroelastic_S22.txt');
time_consol_air_abaqus=(1e3)*sim(:,1); %convert to ms
S22_consol_air_abaqus=scale_stress*(1e-6)*sim(:,2); %convert to MPa
sim=load('disk_axisymm_abaqus_poroelastic_POR.txt');
POR_consol_air_abaqus=scale_stress*(1e-6)*sim(:,2); %convert to MPa

%%%%%%%%%%%%%%%%%%%%%%%%%read nodal tahoe data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extract variable pore_press (y/n) ? y
%  Extract variable d_x (y/n) ? n
%  Extract variable d_y (y/n) ? n
%  Extract variable d_z (y/n) ? y
%  Extract variable s11 (y/n) ? n
%  Extract variable s22 (y/n) ? n
%  Extract variable s33 (y/n) ? y
%  Extract variable s23 (y/n) ? n
%  Extract variable s13 (y/n) ? n
%  Extract variable s12 (y/n) ? n
%  Extract variable p_f (y/n) ? y
%  Extract variable e11 (y/n) ? n
%  Extract variable e22 (y/n) ? n
%  Extract variable e33 (y/n) ? y
%  Extract variable e23 (y/n) ? n
%  Extract variable e13 (y/n) ? n
%  Extract variable e12 (y/n) ? n
%  Extract variable phi_s (y/n) ? n
%  Extract variable phi_f (y/n) ? y
%  Extract variable J (y/n) ? y
%  Extract variable k (y/n) ? n
%  Extract variable kappa (y/n) ? n
%  Extract variable c (y/n) ? n
%  Extract variable p_prime (y/n) ? y
%  Extract variable sdev_sdev (y/n) ? n
%  Extract variable eps_vol_p (y/n) ? n

%%%%%%%%%%%%pore air%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('nddata-QuarterCyl-Drained-Displ_76.txt');
time_drained_air_76=sim(:,1);
pore_press_drained_air_76=scale_stress*sim(:,2);
dz_drained_air_76=sim(:,3);
s33_drained_air_76=scale_stress*sim(:,4);
pf_drained_air_76=scale_stress*sim(:,5);
e33_drained_air_76=sim(:,6);
phi_f_drained_air_76=sim(:,7);
J_drained_air_76=sim(:,8);
p_prime_drained_air_76=scale_stress*sim(:,9);

sim=load('nddata-QuarterCyl-Consol-Displ_76.txt');
time_consol_air_76=sim(:,1);
pore_press_consol_air_76=scale_stress*sim(:,2);
dz_consol_air_76=sim(:,3);
s33_consol_air_76=scale_stress*sim(:,4);
pf_consol_air_76=scale_stress*sim(:,5);
e33_consol_air_76=sim(:,6);
phi_f_consol_air_76=sim(:,7);
J_consol_air_76=sim(:,8);
p_prime_consol_air_76=scale_stress*sim(:,9);

sim=load('nddata-QuarterCyl-Dyn-Displ_76.txt');
time_dyn_air_76=sim(:,1);
pore_press_dyn_air_76=scale_stress*sim(:,2);
dz_dyn_air_76=sim(:,3);
s33_dyn_air_76=scale_stress*sim(:,4);
pf_dyn_air_76=scale_stress*sim(:,5);
e33_dyn_air_76=sim(:,6);
phi_f_dyn_air_76=sim(:,7);
J_dyn_air_76=sim(:,8);
p_prime_dyn_air_76=scale_stress*sim(:,9);

% stress, sigmazz
figure(3)
plot(time_drained_air_76,s33_drained_air_76,'-k',...
     time_drained_air_abaqus,S22_drained_air_abaqus ,'--k',...
     time_consol_air_76,s33_consol_air_76,'-b',...
     time_consol_air_abaqus,S22_consol_air_abaqus,'--b',...
     time_dyn_air_76,s33_dyn_air_76,'-r',...
    'LineWidth',1.5)
xlabel('TIME (ms)')
ylabel('STRESS (kPa)')
grid on
title('center node 76 - pore air')
legend('\sigma_{zz}^\prime drained','S22 drained - Abaqus',...
       '\sigma_{zz}^\prime consol','S22 consol - Abaqus',...
       '\sigma_{zz}^\prime dyn'...
      )
set(gca,'FontName','Helvetica','FontSize',16)

% stress, pf
figure(4)
plot(time_drained_air_76,pf_drained_air_76,'o-k',...
     time_consol_air_76,pf_consol_air_76,'o-b',...
     time_consol_air_abaqus,POR_consol_air_abaqus,'o--b',...
     time_dyn_air_76,pf_dyn_air_76,'o-r',...
    'LineWidth',1.5)
xlabel('TIME (ms)')
ylabel('STRESS (kPa)')
grid on
title('center node 76 - pore air')
legend('p_f drained',...
       'p_f consol','POR consol - Abaqus',...
       'p_f dyn'...
      )
set(gca,'FontName','Helvetica','FontSize',16)

% % stress, strain
% figure(4)
% plot(-e33_drained_air_91,-s33_drained_air_91,'-k',...
%      -e33_consol_air_91,-s33_consol_air_91,'-b',...
%      -e33_impulse_air_91,-s33_impulse_air_91,'-r',...
%     'LineWidth',1.5)
% xlabel('-AXIAL STRAIN (m/m)')
% ylabel('-AXIAL STRESS (Pa)')
% grid on
% title('mid-length corner node 91 - pore air')
% legend('\sigma_{zz} drained','\sigma_{zz} consol','\sigma_{zz} impulse')
% set(gca,'FontName','Helvetica','FontSize',16)

