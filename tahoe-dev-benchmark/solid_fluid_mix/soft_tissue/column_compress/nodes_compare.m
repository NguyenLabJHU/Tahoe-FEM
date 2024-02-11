
%close all
clear all

format short e

%  Extract variable pore_press (y/n) ? y
%  Extract variable d_x (y/n) ? y
%  Extract variable d_y (y/n) ? y
%  Extract variable d_z (y/n) ? y
%  Extract variable s11 (y/n) ? y
%  Extract variable s22 (y/n) ? y
%  Extract variable s33 (y/n) ? y
%  Extract variable s23 (y/n) ? y
%  Extract variable s13 (y/n) ? y
%  Extract variable s12 (y/n) ? y
%  Extract variable p_f (y/n) ? y
%  Extract variable e11 (y/n) ? y
%  Extract variable e22 (y/n) ? y
%  Extract variable e33 (y/n) ? y
%  Extract variable e23 (y/n) ? y
%  Extract variable e13 (y/n) ? y
%  Extract variable e12 (y/n) ? y
%  Extract variable phi_s (y/n) ? y
%  Extract variable phi_f (y/n) ? y
%  Extract variable J (y/n) ? y
%  Extract variable k (y/n) ? y
%  Extract variable kappa (y/n) ? n
%  Extract variable c (y/n) ? n
%  Extract variable p_prime (y/n) ? y
%  Extract variable sdev_sdev (y/n) ? y
%  Extract variable eps_vol_p (y/n) ? n

%%%%%%%%%%%%pore air%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('nddata_drained_air_08.txt');
time_drained_air_08=sim(:,1);
pore_press_drained_air_08=sim(:,2);
dx_drained_air_08=sim(:,3);
dy_drained_air_08=sim(:,4);
dz_drained_air_08=sim(:,5);
s11_drained_air_08=sim(:,6);
s22_drained_air_08=sim(:,7);
s33_drained_air_08=sim(:,8);
s23_drained_air_08=sim(:,9);
s13_drained_air_08=sim(:,10);
s12_drained_air_08=sim(:,11);
pf_drained_air_08=sim(:,12);
e11_drained_air_08=sim(:,13);
e22_drained_air_08=sim(:,14);
e33_drained_air_08=sim(:,15);
e23_drained_air_08=sim(:,16);
e13_drained_air_08=sim(:,17);
e12_drained_air_08=sim(:,18);
phi_s_drained_air_08=sim(:,19);
phi_f_drained_air_08=sim(:,20);
J_drained_air_08=sim(:,21);
k_drained_air_08=sim(:,22);
p_prime_drained_air_08=sim(:,23);
sdev_sdev_drained_air_08=sim(:,24);

sim=load('nddata_drained_air_91.txt');
time_drained_air_91=sim(:,1);
pore_press_drained_air_91=sim(:,2);
dx_drained_air_91=sim(:,3);
dy_drained_air_91=sim(:,4);
dz_drained_air_91=sim(:,5);
s11_drained_air_91=sim(:,6);
s22_drained_air_91=sim(:,7);
s33_drained_air_91=sim(:,8);
s23_drained_air_91=sim(:,9);
s13_drained_air_91=sim(:,10);
s12_drained_air_91=sim(:,11);
pf_drained_air_91=sim(:,12);
e11_drained_air_91=sim(:,13);
e22_drained_air_91=sim(:,14);
e33_drained_air_91=sim(:,15);
e23_drained_air_91=sim(:,16);
e13_drained_air_91=sim(:,17);
e12_drained_air_91=sim(:,18);
phi_s_drained_air_91=sim(:,19);
phi_f_drained_air_91=sim(:,20);
J_drained_air_91=sim(:,21);
k_drained_air_91=sim(:,22);
p_prime_drained_air_91=sim(:,23);
sdev_sdev_drained_air_91=sim(:,24);

sim=load('nddata_consol_air_08.txt');
time_consol_air_08=sim(:,1);
pore_press_consol_air_08=sim(:,2);
dx_consol_air_08=sim(:,3);
dy_consol_air_08=sim(:,4);
dz_consol_air_08=sim(:,5);
s11_consol_air_08=sim(:,6);
s22_consol_air_08=sim(:,7);
s33_consol_air_08=sim(:,8);
s23_consol_air_08=sim(:,9);
s13_consol_air_08=sim(:,10);
s12_consol_air_08=sim(:,11);
pf_consol_air_08=sim(:,12);
e11_consol_air_08=sim(:,13);
e22_consol_air_08=sim(:,14);
e33_consol_air_08=sim(:,15);
e23_consol_air_08=sim(:,16);
e13_consol_air_08=sim(:,17);
e12_consol_air_08=sim(:,18);
phi_s_consol_air_08=sim(:,19);
phi_f_consol_air_08=sim(:,20);
J_consol_air_08=sim(:,21);
k_consol_air_08=sim(:,22);
p_prime_consol_air_08=sim(:,23);
sdev_sdev_consol_air_08=sim(:,24);

sim=load('nddata_consol_air_91.txt');
time_consol_air_91=sim(:,1);
pore_press_consol_air_91=sim(:,2);
dx_consol_air_91=sim(:,3);
dy_consol_air_91=sim(:,4);
dz_consol_air_91=sim(:,5);
s11_consol_air_91=sim(:,6);
s22_consol_air_91=sim(:,7);
s33_consol_air_91=sim(:,8);
s23_consol_air_91=sim(:,9);
s13_consol_air_91=sim(:,10);
s12_consol_air_91=sim(:,11);
pf_consol_air_91=sim(:,12);
e11_consol_air_91=sim(:,13);
e22_consol_air_91=sim(:,14);
e33_consol_air_91=sim(:,15);
e23_consol_air_91=sim(:,16);
e13_consol_air_91=sim(:,17);
e12_consol_air_91=sim(:,18);
phi_s_consol_air_91=sim(:,19);
phi_f_consol_air_91=sim(:,20);
J_consol_air_91=sim(:,21);
k_consol_air_91=sim(:,22);
p_prime_consol_air_91=sim(:,23);
sdev_sdev_consol_air_91=sim(:,24);

sim=load('nddata_impulse_air_08.txt');
time_impulse_air_08=sim(:,1);
pore_press_impulse_air_08=sim(:,2);
dx_impulse_air_08=sim(:,3);
dy_impulse_air_08=sim(:,4);
dz_impulse_air_08=sim(:,5);
s11_impulse_air_08=sim(:,6);
s22_impulse_air_08=sim(:,7);
s33_impulse_air_08=sim(:,8);
s23_impulse_air_08=sim(:,9);
s13_impulse_air_08=sim(:,10);
s12_impulse_air_08=sim(:,11);
pf_impulse_air_08=sim(:,12);
e11_impulse_air_08=sim(:,13);
e22_impulse_air_08=sim(:,14);
e33_impulse_air_08=sim(:,15);
e23_impulse_air_08=sim(:,16);
e13_impulse_air_08=sim(:,17);
e12_impulse_air_08=sim(:,18);
phi_s_impulse_air_08=sim(:,19);
phi_f_impulse_air_08=sim(:,20);
J_impulse_air_08=sim(:,21);
k_impulse_air_08=sim(:,22);
p_prime_impulse_air_08=sim(:,23);
sdev_sdev_impulse_air_08=sim(:,24);

sim=load('nddata_impulse_air_91.txt');
time_impulse_air_91=sim(:,1);
pore_press_impulse_air_91=sim(:,2);
dx_impulse_air_91=sim(:,3);
dy_impulse_air_91=sim(:,4);
dz_impulse_air_91=sim(:,5);
s11_impulse_air_91=sim(:,6);
s22_impulse_air_91=sim(:,7);
s33_impulse_air_91=sim(:,8);
s23_impulse_air_91=sim(:,9);
s13_impulse_air_91=sim(:,10);
s12_impulse_air_91=sim(:,11);
pf_impulse_air_91=sim(:,12);
e11_impulse_air_91=sim(:,13);
e22_impulse_air_91=sim(:,14);
e33_impulse_air_91=sim(:,15);
e23_impulse_air_91=sim(:,16);
e13_impulse_air_91=sim(:,17);
e12_impulse_air_91=sim(:,18);
phi_s_impulse_air_91=sim(:,19);
phi_f_impulse_air_91=sim(:,20);
J_impulse_air_91=sim(:,21);
k_impulse_air_91=sim(:,22);
p_prime_impulse_air_91=sim(:,23);
sdev_sdev_impulse_air_91=sim(:,24);

sim=load('nddata_impulse_air_undrained_08.txt');
time_impulse_air_undrained_08=sim(:,1);
pore_press_impulse_air_undrained_08=sim(:,2);
dx_impulse_air_undrained_08=sim(:,3);
dy_impulse_air_undrained_08=sim(:,4);
dz_impulse_air_undrained_08=sim(:,5);
s11_impulse_air_undrained_08=sim(:,6);
s22_impulse_air_undrained_08=sim(:,7);
s33_impulse_air_undrained_08=sim(:,8);
s23_impulse_air_undrained_08=sim(:,9);
s13_impulse_air_undrained_08=sim(:,10);
s12_impulse_air_undrained_08=sim(:,11);
pf_impulse_air_undrained_08=sim(:,12);
e11_impulse_air_undrained_08=sim(:,13);
e22_impulse_air_undrained_08=sim(:,14);
e33_impulse_air_undrained_08=sim(:,15);
e23_impulse_air_undrained_08=sim(:,16);
e13_impulse_air_undrained_08=sim(:,17);
e12_impulse_air_undrained_08=sim(:,18);
phi_s_impulse_air_undrained_08=sim(:,19);
phi_f_impulse_air_undrained_08=sim(:,20);
J_impulse_air_undrained_08=sim(:,21);
k_impulse_air_undrained_08=sim(:,22);
p_prime_impulse_air_undrained_08=sim(:,23);
sdev_sdev_impulse_air_undrained_08=sim(:,24);

sim=load('nddata_impulse_air_undrained_91.txt');
time_impulse_air_undrained_91=sim(:,1);
pore_press_impulse_air_undrained_91=sim(:,2);
dx_impulse_air_undrained_91=sim(:,3);
dy_impulse_air_undrained_91=sim(:,4);
dz_impulse_air_undrained_91=sim(:,5);
s11_impulse_air_undrained_91=sim(:,6);
s22_impulse_air_undrained_91=sim(:,7);
s33_impulse_air_undrained_91=sim(:,8);
s23_impulse_air_undrained_91=sim(:,9);
s13_impulse_air_undrained_91=sim(:,10);
s12_impulse_air_undrained_91=sim(:,11);
pf_impulse_air_undrained_91=sim(:,12);
e11_impulse_air_undrained_91=sim(:,13);
e22_impulse_air_undrained_91=sim(:,14);
e33_impulse_air_undrained_91=sim(:,15);
e23_impulse_air_undrained_91=sim(:,16);
e13_impulse_air_undrained_91=sim(:,17);
e12_impulse_air_undrained_91=sim(:,18);
phi_s_impulse_air_undrained_91=sim(:,19);
phi_f_impulse_air_undrained_91=sim(:,20);
J_impulse_air_undrained_91=sim(:,21);
k_impulse_air_undrained_91=sim(:,22);
p_prime_impulse_air_undrained_91=sim(:,23);
sdev_sdev_impulse_air_undrained_91=sim(:,24);

% % node displ time
% figure(1)
% plot(time_drained_air_08,dz_drained_air_08,'-.k',...
%      time_consol_air_08,dz_consol_air_08,'-.b',...
%     'LineWidth',1.5)
%     %time_impulse_air_08,dz_impulse_air_08,'-.r',...
% xlabel('TIME (sec)')
% ylabel('VERTICAL DISPLACEMENT (m)')
% grid on
% title('top corner node 8 - air')
% legend('drained','consolidating')
% %legend('drained','consolidating','dynamic impulse')
% set(gca,'FontName','Helvetica','FontSize',16)
% %set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

% node displ time
figure(11)
plot(time_drained_air_08,dz_drained_air_08,'o-.k',...
     time_consol_air_08,dz_consol_air_08,'o-.b',...
     time_drained_air_91,dz_drained_air_91,'x-.k',...
     time_consol_air_91,dz_consol_air_91,'x-.b',...
    'LineWidth',1.5)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
grid on
title('pore air')
legend('drained - nd8','consolidating - nd8','drained - nd91','consolidating - nd91')
set(gca,'FontName','Helvetica','FontSize',16)

% % node displ time
% figure(2)
% plot(time_drained_air_91,dz_drained_air_91,'-.k',...
%      time_consol_air_91,dz_consol_air_91,'-.b',...
%     'LineWidth',1.5)
%     %time_impulse_air_91,dz_impulse_air_91,'-.r',...
% xlabel('TIME (sec)')
% ylabel('VERTICAL DISPLACEMENT (m)')
% grid on
% title('mid-length corner node 91 - air')
% legend('drained','consolidating')
% %legend('drained','consolidating','dynamic impulse')
% set(gca,'FontName','Helvetica','FontSize',16)

% node displ time
figure(22)
plot(time_impulse_air_08,dz_impulse_air_08,'o-.k',...
     time_impulse_air_91,dz_impulse_air_91,'x--k',...
     time_impulse_air_undrained_08,dz_impulse_air_undrained_08,'o-.r',...
     time_impulse_air_undrained_91,dz_impulse_air_undrained_91,'x--r',...
    'LineWidth',1.5)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
grid on
title('dynamic impulse - pore air')
legend('node 8','node 91','node 8 - undrained','node 91 - undrained')
set(gca,'FontName','Helvetica','FontSize',16)

% stress, pf
figure(3)
plot(time_drained_air_91,pf_drained_air_91,'o-.k',...
     time_drained_air_91,s33_drained_air_91,'-k',...
     time_drained_air_91,p_prime_drained_air_91,'x--k',...
     time_consol_air_91,pf_consol_air_91,'o-.b',...
     time_consol_air_91,s33_consol_air_91,'-b',...
     time_consol_air_91,p_prime_consol_air_91,'x--b',...
    'LineWidth',1.5)
%      time_impulse_air_91,pf_impulse_air_91,'-.r',...
%      time_impulse_air_91,s33_impulse_air_91,'-r',...
%      time_impulse_air_91,p_prime_impulse_air_91,'--r',...
xlabel('TIME (sec)')
ylabel('STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore air')
legend('p_f drained','\sigma_{zz} drained','p^\prime drained',...
       'p_f consol','\sigma_{zz} consol','p^\prime consol')
       %'p_f impulse','\sigma_{zz} impulse','p^\prime impulse')
set(gca,'FontName','Helvetica','FontSize',16)

% stress, pf
figure(33)
plot(time_impulse_air_91,pf_impulse_air_91,'o-.k',...
     time_impulse_air_91,s33_impulse_air_91,'x-k',...
     time_impulse_air_91,p_prime_impulse_air_91,'--k',...
     time_impulse_air_undrained_91,pf_impulse_air_undrained_91,'o-.r',...
     time_impulse_air_undrained_91,s33_impulse_air_undrained_91,'x-r',...
     time_impulse_air_undrained_91,p_prime_impulse_air_undrained_91,'--r',...
    'LineWidth',1.5)
xlabel('TIME (sec)')
ylabel('STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore air')
legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
       'p_f impulse - undr','\sigma_{zz} impulse - undr','p^\prime impulse - undr')
set(gca,'FontName','Helvetica','FontSize',16)

% stress, strain
figure(4)
plot(-e33_drained_air_91,-s33_drained_air_91,'-k',...
     -e33_consol_air_91,-s33_consol_air_91,'-b',...
     -e33_impulse_air_91,-s33_impulse_air_91,'-r',...
    'LineWidth',1.5)
xlabel('-AXIAL STRAIN (m/m)')
ylabel('-AXIAL STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore air')
legend('\sigma_{zz} drained','\sigma_{zz} consol','\sigma_{zz} impulse')
set(gca,'FontName','Helvetica','FontSize',16)

%%%%%%%%%%%%pore water%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('nddata_drained_water_08.txt');
time_drained_water_08=sim(:,1);
pore_press_drained_water_08=sim(:,2);
dx_drained_water_08=sim(:,3);
dy_drained_water_08=sim(:,4);
dz_drained_water_08=sim(:,5);
s11_drained_water_08=sim(:,6);
s22_drained_water_08=sim(:,7);
s33_drained_water_08=sim(:,8);
s23_drained_water_08=sim(:,9);
s13_drained_water_08=sim(:,10);
s12_drained_water_08=sim(:,11);
pf_drained_water_08=sim(:,12);
e11_drained_water_08=sim(:,13);
e22_drained_water_08=sim(:,14);
e33_drained_water_08=sim(:,15);
e23_drained_water_08=sim(:,16);
e13_drained_water_08=sim(:,17);
e12_drained_water_08=sim(:,18);
phi_s_drained_water_08=sim(:,19);
phi_f_drained_water_08=sim(:,20);
J_drained_water_08=sim(:,21);
k_drained_water_08=sim(:,22);
p_prime_drained_water_08=sim(:,23);
sdev_sdev_drained_water_08=sim(:,24);

sim=load('nddata_drained_water_91.txt');
time_drained_water_91=sim(:,1);
pore_press_drained_water_91=sim(:,2);
dx_drained_water_91=sim(:,3);
dy_drained_water_91=sim(:,4);
dz_drained_water_91=sim(:,5);
s11_drained_water_91=sim(:,6);
s22_drained_water_91=sim(:,7);
s33_drained_water_91=sim(:,8);
s23_drained_water_91=sim(:,9);
s13_drained_water_91=sim(:,10);
s12_drained_water_91=sim(:,11);
pf_drained_water_91=sim(:,12);
e11_drained_water_91=sim(:,13);
e22_drained_water_91=sim(:,14);
e33_drained_water_91=sim(:,15);
e23_drained_water_91=sim(:,16);
e13_drained_water_91=sim(:,17);
e12_drained_water_91=sim(:,18);
phi_s_drained_water_91=sim(:,19);
phi_f_drained_water_91=sim(:,20);
J_drained_water_91=sim(:,21);
k_drained_water_91=sim(:,22);
p_prime_drained_water_91=sim(:,23);
sdev_sdev_drained_water_91=sim(:,24);

sim=load('nddata_drained_water_natBC_08.txt');
time_drained_water_natBC_08=sim(:,1);
pore_press_drained_water_natBC_08=sim(:,2);
dx_drained_water_natBC_08=sim(:,3);
dy_drained_water_natBC_08=sim(:,4);
dz_drained_water_natBC_08=sim(:,5);
s11_drained_water_natBC_08=sim(:,6);
s22_drained_water_natBC_08=sim(:,7);
s33_drained_water_natBC_08=sim(:,8);
s23_drained_water_natBC_08=sim(:,9);
s13_drained_water_natBC_08=sim(:,10);
s12_drained_water_natBC_08=sim(:,11);
pf_drained_water_natBC_08=sim(:,12);
e11_drained_water_natBC_08=sim(:,13);
e22_drained_water_natBC_08=sim(:,14);
e33_drained_water_natBC_08=sim(:,15);
e23_drained_water_natBC_08=sim(:,16);
e13_drained_water_natBC_08=sim(:,17);
e12_drained_water_natBC_08=sim(:,18);
phi_s_drained_water_natBC_08=sim(:,19);
phi_f_drained_water_natBC_08=sim(:,20);
J_drained_water_natBC_08=sim(:,21);
k_drained_water_natBC_08=sim(:,22);
p_prime_drained_water_natBC_08=sim(:,23);
sdev_sdev_drained_water_natBC_08=sim(:,24);

sim=load('nddata_drained_water_natBC_91.txt');
time_drained_water_natBC_91=sim(:,1);
pore_press_drained_water_natBC_91=sim(:,2);
dx_drained_water_natBC_91=sim(:,3);
dy_drained_water_natBC_91=sim(:,4);
dz_drained_water_natBC_91=sim(:,5);
s11_drained_water_natBC_91=sim(:,6);
s22_drained_water_natBC_91=sim(:,7);
s33_drained_water_natBC_91=sim(:,8);
s23_drained_water_natBC_91=sim(:,9);
s13_drained_water_natBC_91=sim(:,10);
s12_drained_water_natBC_91=sim(:,11);
pf_drained_water_natBC_91=sim(:,12);
e11_drained_water_natBC_91=sim(:,13);
e22_drained_water_natBC_91=sim(:,14);
e33_drained_water_natBC_91=sim(:,15);
e23_drained_water_natBC_91=sim(:,16);
e13_drained_water_natBC_91=sim(:,17);
e12_drained_water_natBC_91=sim(:,18);
phi_s_drained_water_natBC_91=sim(:,19);
phi_f_drained_water_natBC_91=sim(:,20);
J_drained_water_natBC_91=sim(:,21);
k_drained_water_natBC_91=sim(:,22);
p_prime_drained_water_natBC_91=sim(:,23);
sdev_sdev_drained_water_natBC_91=sim(:,24);

sim=load('nddata_consol_water_08.txt');
time_consol_water_08=sim(:,1);
pore_press_consol_water_08=sim(:,2);
dx_consol_water_08=sim(:,3);
dy_consol_water_08=sim(:,4);
dz_consol_water_08=sim(:,5);
s11_consol_water_08=sim(:,6);
s22_consol_water_08=sim(:,7);
s33_consol_water_08=sim(:,8);
s23_consol_water_08=sim(:,9);
s13_consol_water_08=sim(:,10);
s12_consol_water_08=sim(:,11);
pf_consol_water_08=sim(:,12);
e11_consol_water_08=sim(:,13);
e22_consol_water_08=sim(:,14);
e33_consol_water_08=sim(:,15);
e23_consol_water_08=sim(:,16);
e13_consol_water_08=sim(:,17);
e12_consol_water_08=sim(:,18);
phi_s_consol_water_08=sim(:,19);
phi_f_consol_water_08=sim(:,20);
J_consol_water_08=sim(:,21);
k_consol_water_08=sim(:,22);
p_prime_consol_water_08=sim(:,23);
sdev_sdev_consol_water_08=sim(:,24);

sim=load('nddata_consol_water_91.txt');
time_consol_water_91=sim(:,1);
pore_press_consol_water_91=sim(:,2);
dx_consol_water_91=sim(:,3);
dy_consol_water_91=sim(:,4);
dz_consol_water_91=sim(:,5);
s11_consol_water_91=sim(:,6);
s22_consol_water_91=sim(:,7);
s33_consol_water_91=sim(:,8);
s23_consol_water_91=sim(:,9);
s13_consol_water_91=sim(:,10);
s12_consol_water_91=sim(:,11);
pf_consol_water_91=sim(:,12);
e11_consol_water_91=sim(:,13);
e22_consol_water_91=sim(:,14);
e33_consol_water_91=sim(:,15);
e23_consol_water_91=sim(:,16);
e13_consol_water_91=sim(:,17);
e12_consol_water_91=sim(:,18);
phi_s_consol_water_91=sim(:,19);
phi_f_consol_water_91=sim(:,20);
J_consol_water_91=sim(:,21);
k_consol_water_91=sim(:,22);
p_prime_consol_water_91=sim(:,23);
sdev_sdev_consol_water_91=sim(:,24);

sim=load('nddata_impulse_water_08.txt');
time_impulse_water_08=sim(:,1);
pore_press_impulse_water_08=sim(:,2);
dx_impulse_water_08=sim(:,3);
dy_impulse_water_08=sim(:,4);
dz_impulse_water_08=sim(:,5);
s11_impulse_water_08=sim(:,6);
s22_impulse_water_08=sim(:,7);
s33_impulse_water_08=sim(:,8);
s23_impulse_water_08=sim(:,9);
s13_impulse_water_08=sim(:,10);
s12_impulse_water_08=sim(:,11);
pf_impulse_water_08=sim(:,12);
e11_impulse_water_08=sim(:,13);
e22_impulse_water_08=sim(:,14);
e33_impulse_water_08=sim(:,15);
e23_impulse_water_08=sim(:,16);
e13_impulse_water_08=sim(:,17);
e12_impulse_water_08=sim(:,18);
phi_s_impulse_water_08=sim(:,19);
phi_f_impulse_water_08=sim(:,20);
J_impulse_water_08=sim(:,21);
k_impulse_water_08=sim(:,22);
p_prime_impulse_water_08=sim(:,23);
sdev_sdev_impulse_water_08=sim(:,24);

sim=load('nddata_impulse_water_91.txt');
time_impulse_water_91=sim(:,1);
pore_press_impulse_water_91=sim(:,2);
dx_impulse_water_91=sim(:,3);
dy_impulse_water_91=sim(:,4);
dz_impulse_water_91=sim(:,5);
s11_impulse_water_91=sim(:,6);
s22_impulse_water_91=sim(:,7);
s33_impulse_water_91=sim(:,8);
s23_impulse_water_91=sim(:,9);
s13_impulse_water_91=sim(:,10);
s12_impulse_water_91=sim(:,11);
pf_impulse_water_91=sim(:,12);
e11_impulse_water_91=sim(:,13);
e22_impulse_water_91=sim(:,14);
e33_impulse_water_91=sim(:,15);
e23_impulse_water_91=sim(:,16);
e13_impulse_water_91=sim(:,17);
e12_impulse_water_91=sim(:,18);
phi_s_impulse_water_91=sim(:,19);
phi_f_impulse_water_91=sim(:,20);
J_impulse_water_91=sim(:,21);
k_impulse_water_91=sim(:,22);
p_prime_impulse_water_91=sim(:,23);
sdev_sdev_impulse_water_91=sim(:,24);

sim=load('nddata_impulse_water_short_08.txt');
time_impulse_water_short_08=sim(:,1);
pore_press_impulse_water_short_08=sim(:,2);
dx_impulse_water_short_08=sim(:,3);
dy_impulse_water_short_08=sim(:,4);
dz_impulse_water_short_08=sim(:,5);
s11_impulse_water_short_08=sim(:,6);
s22_impulse_water_short_08=sim(:,7);
s33_impulse_water_short_08=sim(:,8);
s23_impulse_water_short_08=sim(:,9);
s13_impulse_water_short_08=sim(:,10);
s12_impulse_water_short_08=sim(:,11);
pf_impulse_water_short_08=sim(:,12);
e11_impulse_water_short_08=sim(:,13);
e22_impulse_water_short_08=sim(:,14);
e33_impulse_water_short_08=sim(:,15);
e23_impulse_water_short_08=sim(:,16);
e13_impulse_water_short_08=sim(:,17);
e12_impulse_water_short_08=sim(:,18);
phi_s_impulse_water_short_08=sim(:,19);
phi_f_impulse_water_short_08=sim(:,20);
J_impulse_water_short_08=sim(:,21);
k_impulse_water_short_08=sim(:,22);
p_prime_impulse_water_short_08=sim(:,23);
sdev_sdev_impulse_water_short_08=sim(:,24);

sim=load('nddata_impulse_water_short_91.txt');
time_impulse_water_short_91=sim(:,1);
pore_press_impulse_water_short_91=sim(:,2);
dx_impulse_water_short_91=sim(:,3);
dy_impulse_water_short_91=sim(:,4);
dz_impulse_water_short_91=sim(:,5);
s11_impulse_water_short_91=sim(:,6);
s22_impulse_water_short_91=sim(:,7);
s33_impulse_water_short_91=sim(:,8);
s23_impulse_water_short_91=sim(:,9);
s13_impulse_water_short_91=sim(:,10);
s12_impulse_water_short_91=sim(:,11);
pf_impulse_water_short_91=sim(:,12);
e11_impulse_water_short_91=sim(:,13);
e22_impulse_water_short_91=sim(:,14);
e33_impulse_water_short_91=sim(:,15);
e23_impulse_water_short_91=sim(:,16);
e13_impulse_water_short_91=sim(:,17);
e12_impulse_water_short_91=sim(:,18);
phi_s_impulse_water_short_91=sim(:,19);
phi_f_impulse_water_short_91=sim(:,20);
J_impulse_water_short_91=sim(:,21);
k_impulse_water_short_91=sim(:,22);
p_prime_impulse_water_short_91=sim(:,23);
sdev_sdev_impulse_water_short_91=sim(:,24);

sim=load('nddata_impulse_water_d_08.txt');
time_impulse_water_d_08=sim(:,1);
pore_press_impulse_water_d_08=sim(:,2);
dx_impulse_water_d_08=sim(:,3);
dy_impulse_water_d_08=sim(:,4);
dz_impulse_water_d_08=sim(:,5);
s11_impulse_water_d_08=sim(:,6);
s22_impulse_water_d_08=sim(:,7);
s33_impulse_water_d_08=sim(:,8);
s23_impulse_water_d_08=sim(:,9);
s13_impulse_water_d_08=sim(:,10);
s12_impulse_water_d_08=sim(:,11);
pf_impulse_water_d_08=sim(:,12);
e11_impulse_water_d_08=sim(:,13);
e22_impulse_water_d_08=sim(:,14);
e33_impulse_water_d_08=sim(:,15);
e23_impulse_water_d_08=sim(:,16);
e13_impulse_water_d_08=sim(:,17);
e12_impulse_water_d_08=sim(:,18);
phi_s_impulse_water_d_08=sim(:,19);
phi_f_impulse_water_d_08=sim(:,20);
J_impulse_water_d_08=sim(:,21);
k_impulse_water_d_08=sim(:,22);
p_prime_impulse_water_d_08=sim(:,23);
sdev_sdev_impulse_water_d_08=sim(:,24);

sim=load('nddata_impulse_water_d_91.txt');
time_impulse_water_d_91=sim(:,1);
pore_press_impulse_water_d_91=sim(:,2);
dx_impulse_water_d_91=sim(:,3);
dy_impulse_water_d_91=sim(:,4);
dz_impulse_water_d_91=sim(:,5);
s11_impulse_water_d_91=sim(:,6);
s22_impulse_water_d_91=sim(:,7);
s33_impulse_water_d_91=sim(:,8);
s23_impulse_water_d_91=sim(:,9);
s13_impulse_water_d_91=sim(:,10);
s12_impulse_water_d_91=sim(:,11);
pf_impulse_water_d_91=sim(:,12);
e11_impulse_water_d_91=sim(:,13);
e22_impulse_water_d_91=sim(:,14);
e33_impulse_water_d_91=sim(:,15);
e23_impulse_water_d_91=sim(:,16);
e13_impulse_water_d_91=sim(:,17);
e12_impulse_water_d_91=sim(:,18);
phi_s_impulse_water_d_91=sim(:,19);
phi_f_impulse_water_d_91=sim(:,20);
J_impulse_water_d_91=sim(:,21);
k_impulse_water_d_91=sim(:,22);
p_prime_impulse_water_d_91=sim(:,23);
sdev_sdev_impulse_water_d_91=sim(:,24);

sim=load('nddata_impulse_water_undrained_08.txt');
time_impulse_water_undrained_08=sim(:,1);
pore_press_impulse_water_undrained_08=sim(:,2);
dx_impulse_water_undrained_08=sim(:,3);
dy_impulse_water_undrained_08=sim(:,4);
dz_impulse_water_undrained_08=sim(:,5);
s11_impulse_water_undrained_08=sim(:,6);
s22_impulse_water_undrained_08=sim(:,7);
s33_impulse_water_undrained_08=sim(:,8);
s23_impulse_water_undrained_08=sim(:,9);
s13_impulse_water_undrained_08=sim(:,10);
s12_impulse_water_undrained_08=sim(:,11);
pf_impulse_water_undrained_08=sim(:,12);
e11_impulse_water_undrained_08=sim(:,13);
e22_impulse_water_undrained_08=sim(:,14);
e33_impulse_water_undrained_08=sim(:,15);
e23_impulse_water_undrained_08=sim(:,16);
e13_impulse_water_undrained_08=sim(:,17);
e12_impulse_water_undrained_08=sim(:,18);
phi_s_impulse_water_undrained_08=sim(:,19);
phi_f_impulse_water_undrained_08=sim(:,20);
J_impulse_water_undrained_08=sim(:,21);
k_impulse_water_undrained_08=sim(:,22);
p_prime_impulse_water_undrained_08=sim(:,23);
sdev_sdev_impulse_water_undrained_08=sim(:,24);

sim=load('nddata_impulse_water_undrained_91.txt');
time_impulse_water_undrained_91=sim(:,1);
pore_press_impulse_water_undrained_91=sim(:,2);
dx_impulse_water_undrained_91=sim(:,3);
dy_impulse_water_undrained_91=sim(:,4);
dz_impulse_water_undrained_91=sim(:,5);
s11_impulse_water_undrained_91=sim(:,6);
s22_impulse_water_undrained_91=sim(:,7);
s33_impulse_water_undrained_91=sim(:,8);
s23_impulse_water_undrained_91=sim(:,9);
s13_impulse_water_undrained_91=sim(:,10);
s12_impulse_water_undrained_91=sim(:,11);
pf_impulse_water_undrained_91=sim(:,12);
e11_impulse_water_undrained_91=sim(:,13);
e22_impulse_water_undrained_91=sim(:,14);
e33_impulse_water_undrained_91=sim(:,15);
e23_impulse_water_undrained_91=sim(:,16);
e13_impulse_water_undrained_91=sim(:,17);
e12_impulse_water_undrained_91=sim(:,18);
phi_s_impulse_water_undrained_91=sim(:,19);
phi_f_impulse_water_undrained_91=sim(:,20);
J_impulse_water_undrained_91=sim(:,21);
k_impulse_water_undrained_91=sim(:,22);
p_prime_impulse_water_undrained_91=sim(:,23);
sdev_sdev_impulse_water_undrained_91=sim(:,24);

sim=load('nddata_impulse_solid_08.txt');
time_impulse_solid_08=sim(:,1);
pore_press_impulse_solid_08=sim(:,2);
dx_impulse_solid_08=sim(:,3);
dy_impulse_solid_08=sim(:,4);
dz_impulse_solid_08=sim(:,5);
s11_impulse_solid_08=sim(:,6);
s22_impulse_solid_08=sim(:,7);
s33_impulse_solid_08=sim(:,8);
s23_impulse_solid_08=sim(:,9);
s13_impulse_solid_08=sim(:,10);
s12_impulse_solid_08=sim(:,11);
pf_impulse_solid_08=sim(:,12);
e11_impulse_solid_08=sim(:,13);
e22_impulse_solid_08=sim(:,14);
e33_impulse_solid_08=sim(:,15);
e23_impulse_solid_08=sim(:,16);
e13_impulse_solid_08=sim(:,17);
e12_impulse_solid_08=sim(:,18);
phi_s_impulse_solid_08=sim(:,19);
phi_f_impulse_solid_08=sim(:,20);
J_impulse_solid_08=sim(:,21);
k_impulse_solid_08=sim(:,22);
p_prime_impulse_solid_08=sim(:,23);
sdev_sdev_impulse_solid_08=sim(:,24);

sim=load('nddata_impulse_solid_91.txt');
time_impulse_solid_91=sim(:,1);
pore_press_impulse_solid_91=sim(:,2);
dx_impulse_solid_91=sim(:,3);
dy_impulse_solid_91=sim(:,4);
dz_impulse_solid_91=sim(:,5);
s11_impulse_solid_91=sim(:,6);
s22_impulse_solid_91=sim(:,7);
s33_impulse_solid_91=sim(:,8);
s23_impulse_solid_91=sim(:,9);
s13_impulse_solid_91=sim(:,10);
s12_impulse_solid_91=sim(:,11);
pf_impulse_solid_91=sim(:,12);
e11_impulse_solid_91=sim(:,13);
e22_impulse_solid_91=sim(:,14);
e33_impulse_solid_91=sim(:,15);
e23_impulse_solid_91=sim(:,16);
e13_impulse_solid_91=sim(:,17);
e12_impulse_solid_91=sim(:,18);
phi_s_impulse_solid_91=sim(:,19);
phi_f_impulse_solid_91=sim(:,20);
J_impulse_solid_91=sim(:,21);
k_impulse_solid_91=sim(:,22);
p_prime_impulse_solid_91=sim(:,23);
sdev_sdev_impulse_solid_91=sim(:,24);

% % node displ time
% figure(5)
% plot(time_drained_water_08,dz_drained_water_08,'-.k',...
%      time_consol_water_08,dz_consol_water_08,'-.b',...
%     'LineWidth',1.5)
%     %time_impulse_water_08,dz_impulse_water_08,'-.r',...
% xlabel('TIME (sec)')
% ylabel('VERTICAL DISPLACEMENT (m)')
% grid on
% title('top corner node 8 - water')
% legend('drained','consolidating')
% %legend('drained','consolidating','dynamic impulse')
% set(gca,'FontName','Helvetica','FontSize',16)

% node displ time
figure(55)
plot(time_drained_water_08,dz_drained_water_08,'o-.k',...
     time_drained_water_natBC_08,dz_drained_water_natBC_08,'.r',...
     time_consol_water_08,dz_consol_water_08,'o-.b',...
     time_drained_water_91,dz_drained_water_91,'x-.k',...
     time_drained_water_natBC_91,dz_drained_water_natBC_91,'.r',...
     time_consol_water_91,dz_consol_water_91,'x-.b',...
	 'LineWidth',1.5)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
grid on
title('pore water')
legend('drained - nd8','drained, natBC - nd8','consolidating - nd8',...
       'drained - nd91','drained, natBC - nd91','consolidating - nd91')
set(gca,'FontName','Helvetica','FontSize',16)

% % node displ time
% figure(6)
% plot(time_drained_water_91,dz_drained_water_91,'-.k',...
%      time_consol_water_91,dz_consol_water_91,'-.b',...
%     'LineWidth',1.5)
% %     time_impulse_water_91,dz_impulse_water_91,'-.r',...
% %     time_impulse_water_d_91,dz_impulse_water_d_91,'-.g',...
% xlabel('TIME (sec)')
% ylabel('VERTICAL DISPLACEMENT (m)')
% grid on
% title('mid-length corner node 91 - water')
% legend('drained','consolidating')
% %legend('drained','consolidating','dynamic impulse','dynamic impulse d')
% set(gca,'FontName','Helvetica','FontSize',16)

% node displ time
figure(66)
plot(time_impulse_water_08,dz_impulse_water_08,'o-.b',...
     time_impulse_water_91,dz_impulse_water_91,'x--b',...
     time_impulse_water_undrained_08,dz_impulse_water_undrained_08,'o-.r',...
     time_impulse_water_undrained_91,dz_impulse_water_undrained_91,'x--r',...
    'LineWidth',1.5)
%     time_impulse_water_d_91,dz_impulse_water_d_91,'-.g',...
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
grid on
title('dynamic impulse - pore water')
%legend('node 8','node 91')
%legend('node 8','node 91','node 91 d')
legend('node 8','node 91','node 8 - undr','node 91 - undr')
set(gca,'FontName','Helvetica','FontSize',16)

% node displ time
figure(666)
plot(time_impulse_water_08,dz_impulse_water_08,'o-.b',...
     time_impulse_water_91,dz_impulse_water_91,'x--b',...
     time_impulse_solid_08,dz_impulse_solid_08,'o-.k',...
     time_impulse_solid_91,dz_impulse_solid_91,'x--k',...
     time_impulse_water_undrained_08,dz_impulse_water_undrained_08,'o-.r',...
     time_impulse_water_undrained_91,dz_impulse_water_undrained_91,'x--r',...
    'LineWidth',1.5)
%     time_impulse_water_d_91,dz_impulse_water_d_91,'-.g',...
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
grid on
title('dynamic impulse - pore water')
%legend('node 8','node 91')
%legend('node 8','node 91','node 91 d')
legend('node 8','node 91','node 8 solid','node 91 solid','node 8 - undr','node 91 - undr')
set(gca,'FontName','Helvetica','FontSize',16)

% % node displ time
% figure(6666)
% plot(time_impulse_water_08,dz_impulse_water_08,'o-.b',...
%      time_impulse_water_91,dz_impulse_water_91,'x--b',...
%      time_impulse_water_short_08,dz_impulse_water_short_08,'o-.g',...
%      time_impulse_water_short_91,dz_impulse_water_short_91,'x--g',...
%      time_impulse_water_undrained_08,dz_impulse_water_undrained_08,'o-.r',...
%      time_impulse_water_undrained_91,dz_impulse_water_undrained_91,'x--r',...
%     'LineWidth',1.5)
% %     time_impulse_water_d_91,dz_impulse_water_d_91,'-.g',...
% xlabel('TIME (sec)')
% ylabel('VERTICAL DISPLACEMENT (m)')
% grid on
% title('dynamic impulse - water')
% %legend('node 8','node 91')
% %legend('node 8','node 91','node 91 d')
% legend('node 8 - 1e-2','node 91 - 1e-2','node 8 - 1e-3','node 91 - 1e-3',...
%        'node 8 - undr','node 91 - undr')
% set(gca,'FontName','Helvetica','FontSize',16)

% stress, pf
figure(7)
plot(time_drained_water_91,pf_drained_water_91,'o-.k',...
     time_drained_water_91,s33_drained_water_91,'-k',...
     time_drained_water_91,p_prime_drained_water_91,'x--k',...
     time_consol_water_91,pf_consol_water_91,'o-.b',...
     time_consol_water_91,s33_consol_water_91,'-b',...
     time_consol_water_91,p_prime_consol_water_91,'x--b',...
    'LineWidth',1.5)
%     time_impulse_water_91,pf_impulse_water_91,'-.r',...
%      time_impulse_water_91,s33_impulse_water_91,'-r',...
%      time_impulse_water_91,p_prime_impulse_water_91,'--r',...
%      time_impulse_water_d_91,pf_impulse_water_d_91,'-.g',...
%      time_impulse_water_d_91,s33_impulse_water_d_91,'-g',...
%      time_impulse_water_d_91,p_prime_impulse_water_d_91,'--g',...
xlabel('TIME (sec)')
ylabel('STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore water')
legend('p_f drained','\sigma_{zz} drained','p^\prime drained',...
       'p_f consol','\sigma_{zz} consol','p^\prime consol')
% legend('p_f drained','\sigma_{zz} drained','p^\prime drained',...
%        'p_f consol','\sigma_{zz} consol','p^\prime consol',...
%        'p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
%        'p_f impulse d','\sigma_{zz} impulse d','p^\prime impulse d')   
set(gca,'FontName','Helvetica','FontSize',16)

% stress, pf
figure(77)
plot(time_impulse_water_91,pf_impulse_water_91,'o-.b',...
     time_impulse_water_91,s33_impulse_water_91,'-b',...
     time_impulse_water_91,p_prime_impulse_water_91,'x--b',...
     time_impulse_water_undrained_91,pf_impulse_water_undrained_91,'o-.r',...
     time_impulse_water_undrained_91,s33_impulse_water_undrained_91,'-r',...
     time_impulse_water_undrained_91,p_prime_impulse_water_undrained_91,'x--r',...
    'LineWidth',1.5)
%     time_impulse_water_d_91,pf_impulse_water_d_91,'-.g',...
%      time_impulse_water_d_91,s33_impulse_water_d_91,'-g',...
%      time_impulse_water_d_91,p_prime_impulse_water_d_91,'--g',...
xlabel('TIME (sec)')
ylabel('STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore water')
legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse')
% legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
%        'p_f impulse d','\sigma_{zz} impulse d','p^\prime impulse d')
legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
       'p_f impulse - undr','\sigma_{zz} impulse - undr','p^\prime impulse - undr')
set(gca,'FontName','Helvetica','FontSize',16)

% stress, pf
figure(777)
plot(time_impulse_water_91,pf_impulse_water_91,'o-.b',...
     time_impulse_water_91,s33_impulse_water_91,'-b',...
     time_impulse_water_91,p_prime_impulse_water_91,'x--b',...
     time_impulse_solid_91,pf_impulse_solid_91,'o-.k',...
     time_impulse_solid_91,s33_impulse_solid_91,'-k',...
     time_impulse_solid_91,p_prime_impulse_solid_91,'x--k',...
     time_impulse_water_undrained_91,pf_impulse_water_undrained_91,'o-.r',...
     time_impulse_water_undrained_91,s33_impulse_water_undrained_91,'-r',...
     time_impulse_water_undrained_91,p_prime_impulse_water_undrained_91,'x--r',...
    'LineWidth',1.5)
%     time_impulse_water_d_91,pf_impulse_water_d_91,'-.g',...
%      time_impulse_water_d_91,s33_impulse_water_d_91,'-g',...
%      time_impulse_water_d_91,p_prime_impulse_water_d_91,'--g',...
xlabel('TIME (sec)')
ylabel('STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore water')
legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse')
% legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
%        'p_f impulse d','\sigma_{zz} impulse d','p^\prime impulse d')
legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
       'p_f impulse solid','\sigma_{zz} impulse solid','p^\prime impulse solid',...
       'p_f impulse - undr','\sigma_{zz} impulse - undr','p^\prime impulse - undr')
set(gca,'FontName','Helvetica','FontSize',16)

% % stress, pf
% figure(7777)
% plot(time_impulse_water_91,pf_impulse_water_91,'o-.b',...
%      time_impulse_water_91,s33_impulse_water_91,'-b',...
%      time_impulse_water_91,p_prime_impulse_water_91,'x--b',...
%      time_impulse_water_short_91,pf_impulse_water_short_91,'o-.k',...
%      time_impulse_water_short_91,s33_impulse_water_short_91,'-k',...
%      time_impulse_water_short_91,p_prime_impulse_water_short_91,'x--k',...
%      time_impulse_water_undrained_91,pf_impulse_water_undrained_91,'o-.r',...
%      time_impulse_water_undrained_91,s33_impulse_water_undrained_91,'-r',...
%      time_impulse_water_undrained_91,p_prime_impulse_water_undrained_91,'x--r',...
%     'LineWidth',1.5)
% %     time_impulse_water_d_91,pf_impulse_water_d_91,'-.g',...
% %      time_impulse_water_d_91,s33_impulse_water_d_91,'-g',...
% %      time_impulse_water_d_91,p_prime_impulse_water_d_91,'--g',...
% xlabel('TIME (sec)')
% ylabel('STRESS (Pa)')
% grid on
% title('mid-length corner node 91 - water')
% legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse')
% % legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
% %        'p_f impulse d','\sigma_{zz} impulse d','p^\prime impulse d')
% legend('p_f impulse','\sigma_{zz} impulse','p^\prime impulse',...
%        'p_f impulse short','\sigma_{zz} impulse short','p^\prime impulse short',...
%        'p_f impulse - undr','\sigma_{zz} impulse - undr','p^\prime impulse - undr')
% set(gca,'FontName','Helvetica','FontSize',16)



% stress, strain
figure(8)
plot(-e33_drained_water_91,-s33_drained_water_91,'-k',...
     -e33_consol_water_91,-s33_consol_water_91,'-b',...
     -e33_impulse_water_91,-s33_impulse_water_91,'-r',...
     -e33_impulse_water_d_91,-s33_impulse_water_d_91,'-g',...
    'LineWidth',1.5)
xlabel('-AXIAL STRAIN (m/m)')
ylabel('-AXIAL STRESS (Pa)')
grid on
title('mid-length corner node 91 - pore water')
legend('\sigma_{zz} drained','\sigma_{zz} consol','\sigma_{zz} impulse','\sigma_{zz} impulse d')
set(gca,'FontName','Helvetica','FontSize',16)

