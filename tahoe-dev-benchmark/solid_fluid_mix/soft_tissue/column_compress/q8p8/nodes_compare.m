
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

%%%%%%%%%%%%pore water%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('../nddata_consol_water_smallk_08.txt');
time_consol_water_smallk_08=sim(:,1);
pore_press_consol_water_smallk_08=sim(:,2);
dx_consol_water_smallk_08=sim(:,3);
dy_consol_water_smallk_08=sim(:,4);
dz_consol_water_smallk_08=sim(:,5);
s11_consol_water_smallk_08=sim(:,6);
s22_consol_water_smallk_08=sim(:,7);
s33_consol_water_smallk_08=sim(:,8);
s23_consol_water_smallk_08=sim(:,9);
s13_consol_water_smallk_08=sim(:,10);
s12_consol_water_smallk_08=sim(:,11);
pf_consol_water_smallk_08=sim(:,12);
e11_consol_water_smallk_08=sim(:,13);
e22_consol_water_smallk_08=sim(:,14);
e33_consol_water_smallk_08=sim(:,15);
e23_consol_water_smallk_08=sim(:,16);
e13_consol_water_smallk_08=sim(:,17);
e12_consol_water_smallk_08=sim(:,18);
phi_s_consol_water_smallk_08=sim(:,19);
phi_f_consol_water_smallk_08=sim(:,20);
J_consol_water_smallk_08=sim(:,21);
k_consol_water_smallk_08=sim(:,22);
p_prime_consol_water_smallk_08=sim(:,23);
sdev_sdev_consol_water_smallk_08=sim(:,24);

sim=load('../nddata_consol_water_smallk_91.txt');
time_consol_water_smallk_91=sim(:,1);
pore_press_consol_water_smallk_91=sim(:,2);
dx_consol_water_smallk_91=sim(:,3);
dy_consol_water_smallk_91=sim(:,4);
dz_consol_water_smallk_91=sim(:,5);
s11_consol_water_smallk_91=sim(:,6);
s22_consol_water_smallk_91=sim(:,7);
s33_consol_water_smallk_91=sim(:,8);
s23_consol_water_smallk_91=sim(:,9);
s13_consol_water_smallk_91=sim(:,10);
s12_consol_water_smallk_91=sim(:,11);
pf_consol_water_smallk_91=sim(:,12);
e11_consol_water_smallk_91=sim(:,13);
e22_consol_water_smallk_91=sim(:,14);
e33_consol_water_smallk_91=sim(:,15);
e23_consol_water_smallk_91=sim(:,16);
e13_consol_water_smallk_91=sim(:,17);
e12_consol_water_smallk_91=sim(:,18);
phi_s_consol_water_smallk_91=sim(:,19);
phi_f_consol_water_smallk_91=sim(:,20);
J_consol_water_smallk_91=sim(:,21);
k_consol_water_smallk_91=sim(:,22);
p_prime_consol_water_smallk_91=sim(:,23);
sdev_sdev_consol_water_smallk_91=sim(:,24);

%%%%%%%%%%%%pore water - q8p8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('./nddata_consol_water_smallk_q8p8_08.txt');
time_consol_water_smallk_q8p8_08=sim(:,1);
pore_press_consol_water_smallk_q8p8_08=sim(:,2);
dx_consol_water_smallk_q8p8_08=sim(:,3);
dy_consol_water_smallk_q8p8_08=sim(:,4);
dz_consol_water_smallk_q8p8_08=sim(:,5);
s11_consol_water_smallk_q8p8_08=sim(:,6);
s22_consol_water_smallk_q8p8_08=sim(:,7);
s33_consol_water_smallk_q8p8_08=sim(:,8);
s23_consol_water_smallk_q8p8_08=sim(:,9);
s13_consol_water_smallk_q8p8_08=sim(:,10);
s12_consol_water_smallk_q8p8_08=sim(:,11);
pf_consol_water_smallk_q8p8_08=sim(:,12);
e11_consol_water_smallk_q8p8_08=sim(:,13);
e22_consol_water_smallk_q8p8_08=sim(:,14);
e33_consol_water_smallk_q8p8_08=sim(:,15);
e23_consol_water_smallk_q8p8_08=sim(:,16);
e13_consol_water_smallk_q8p8_08=sim(:,17);
e12_consol_water_smallk_q8p8_08=sim(:,18);
phi_s_consol_water_smallk_q8p8_08=sim(:,19);
phi_f_consol_water_smallk_q8p8_08=sim(:,20);
J_consol_water_smallk_q8p8_08=sim(:,21);
k_consol_water_smallk_q8p8_08=sim(:,22);
p_prime_consol_water_smallk_q8p8_08=sim(:,23);
sdev_sdev_consol_water_smallk_q8p8_08=sim(:,24);

sim=load('./nddata_consol_water_smallk_q8p8_24.txt');
time_consol_water_smallk_q8p8_24=sim(:,1);
pore_press_consol_water_smallk_q8p8_24=sim(:,2);
dx_consol_water_smallk_q8p8_24=sim(:,3);
dy_consol_water_smallk_q8p8_24=sim(:,4);
dz_consol_water_smallk_q8p8_24=sim(:,5);
s11_consol_water_smallk_q8p8_24=sim(:,6);
s22_consol_water_smallk_q8p8_24=sim(:,7);
s33_consol_water_smallk_q8p8_24=sim(:,8);
s23_consol_water_smallk_q8p8_24=sim(:,9);
s13_consol_water_smallk_q8p8_24=sim(:,10);
s12_consol_water_smallk_q8p8_24=sim(:,11);
pf_consol_water_smallk_q8p8_24=sim(:,12);
e11_consol_water_smallk_q8p8_24=sim(:,13);
e22_consol_water_smallk_q8p8_24=sim(:,14);
e33_consol_water_smallk_q8p8_24=sim(:,15);
e23_consol_water_smallk_q8p8_24=sim(:,16);
e13_consol_water_smallk_q8p8_24=sim(:,17);
e12_consol_water_smallk_q8p8_24=sim(:,18);
phi_s_consol_water_smallk_q8p8_24=sim(:,19);
phi_f_consol_water_smallk_q8p8_24=sim(:,20);
J_consol_water_smallk_q8p8_24=sim(:,21);
k_consol_water_smallk_q8p8_24=sim(:,22);
p_prime_consol_water_smallk_q8p8_24=sim(:,23);
sdev_sdev_consol_water_smallk_q8p8_24=sim(:,24);

%%%%%%%%%%%%pore water - q8p8 - Alpha1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim=load('./nddata_consol_water_smallk_q8p8_Alpha1_08.txt');
time_consol_water_smallk_q8p8_Alpha1_08=sim(:,1);
pore_press_consol_water_smallk_q8p8_Alpha1_08=sim(:,2);
dx_consol_water_smallk_q8p8_Alpha1_08=sim(:,3);
dy_consol_water_smallk_q8p8_Alpha1_08=sim(:,4);
dz_consol_water_smallk_q8p8_Alpha1_08=sim(:,5);
s11_consol_water_smallk_q8p8_Alpha1_08=sim(:,6);
s22_consol_water_smallk_q8p8_Alpha1_08=sim(:,7);
s33_consol_water_smallk_q8p8_Alpha1_08=sim(:,8);
s23_consol_water_smallk_q8p8_Alpha1_08=sim(:,9);
s13_consol_water_smallk_q8p8_Alpha1_08=sim(:,10);
s12_consol_water_smallk_q8p8_Alpha1_08=sim(:,11);
pf_consol_water_smallk_q8p8_Alpha1_08=sim(:,12);
e11_consol_water_smallk_q8p8_Alpha1_08=sim(:,13);
e22_consol_water_smallk_q8p8_Alpha1_08=sim(:,14);
e33_consol_water_smallk_q8p8_Alpha1_08=sim(:,15);
e23_consol_water_smallk_q8p8_Alpha1_08=sim(:,16);
e13_consol_water_smallk_q8p8_Alpha1_08=sim(:,17);
e12_consol_water_smallk_q8p8_Alpha1_08=sim(:,18);
phi_s_consol_water_smallk_q8p8_Alpha1_08=sim(:,19);
phi_f_consol_water_smallk_q8p8_Alpha1_08=sim(:,20);
J_consol_water_smallk_q8p8_Alpha1_08=sim(:,21);
k_consol_water_smallk_q8p8_Alpha1_08=sim(:,22);
p_prime_consol_water_smallk_q8p8_Alpha1_08=sim(:,23);
sdev_sdev_consol_water_smallk_q8p8_Alpha1_08=sim(:,24);

sim=load('./nddata_consol_water_smallk_q8p8_Alpha1_24.txt');
time_consol_water_smallk_q8p8_Alpha1_24=sim(:,1);
pore_press_consol_water_smallk_q8p8_Alpha1_24=sim(:,2);
dx_consol_water_smallk_q8p8_Alpha1_24=sim(:,3);
dy_consol_water_smallk_q8p8_Alpha1_24=sim(:,4);
dz_consol_water_smallk_q8p8_Alpha1_24=sim(:,5);
s11_consol_water_smallk_q8p8_Alpha1_24=sim(:,6);
s22_consol_water_smallk_q8p8_Alpha1_24=sim(:,7);
s33_consol_water_smallk_q8p8_Alpha1_24=sim(:,8);
s23_consol_water_smallk_q8p8_Alpha1_24=sim(:,9);
s13_consol_water_smallk_q8p8_Alpha1_24=sim(:,10);
s12_consol_water_smallk_q8p8_Alpha1_24=sim(:,11);
pf_consol_water_smallk_q8p8_Alpha1_24=sim(:,12);
e11_consol_water_smallk_q8p8_Alpha1_24=sim(:,13);
e22_consol_water_smallk_q8p8_Alpha1_24=sim(:,14);
e33_consol_water_smallk_q8p8_Alpha1_24=sim(:,15);
e23_consol_water_smallk_q8p8_Alpha1_24=sim(:,16);
e13_consol_water_smallk_q8p8_Alpha1_24=sim(:,17);
e12_consol_water_smallk_q8p8_Alpha1_24=sim(:,18);
phi_s_consol_water_smallk_q8p8_Alpha1_24=sim(:,19);
phi_f_consol_water_smallk_q8p8_Alpha1_24=sim(:,20);
J_consol_water_smallk_q8p8_Alpha1_24=sim(:,21);
k_consol_water_smallk_q8p8_Alpha1_24=sim(:,22);
p_prime_consol_water_smallk_q8p8_Alpha1_24=sim(:,23);
sdev_sdev_consol_water_smallk_q8p8_Alpha1_24=sim(:,24);

% node displ time
figure(55)
plot(time_consol_water_smallk_08,dz_consol_water_smallk_08,'o-.b',...
     time_consol_water_smallk_91,dz_consol_water_smallk_91,'x-.b',...
     time_consol_water_smallk_q8p8_08,dz_consol_water_smallk_q8p8_08,'o-.r',...
     time_consol_water_smallk_q8p8_24,dz_consol_water_smallk_q8p8_24,'x-.r',...
     time_consol_water_smallk_q8p8_Alpha1_08,dz_consol_water_smallk_q8p8_Alpha1_08,'o-.g',...
     time_consol_water_smallk_q8p8_Alpha1_24,dz_consol_water_smallk_q8p8_Alpha1_24,'x-.g',...
     'LineWidth',1.5)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
grid on
title('pore water')
legend('consolidating - nd8',...
       'consolidating - nd91',...
       'consolidating q8p8 - nd8',...
       'consolidating q8p8 - nd24',...
       'consolidating q8p8 Alpha1 - nd8',...
       'consolidating q8p8 Alpha1 - nd24')
set(gca,'FontName','Helvetica','FontSize',16)

% node displ time
figure(66)
plot(time_consol_water_smallk_08,pf_consol_water_smallk_08,'o-.b',...
     time_consol_water_smallk_91,pf_consol_water_smallk_91,'x-.b',...
     time_consol_water_smallk_q8p8_08,pf_consol_water_smallk_q8p8_08,'o-.r',...
     time_consol_water_smallk_q8p8_24,pf_consol_water_smallk_q8p8_24,'x-.r',...
     time_consol_water_smallk_q8p8_Alpha1_08,pf_consol_water_smallk_q8p8_Alpha1_08,'o-.g',...
     time_consol_water_smallk_q8p8_Alpha1_24,pf_consol_water_smallk_q8p8_Alpha1_24,'x-.g',...
     'LineWidth',1.5)
xlabel('TIME (sec)')
ylabel('PORE PRESSURE (Pa)')
grid on
title('pore water')
legend('consolidating - nd8',...
       'consolidating - nd91',...
       'consolidating q8p8 - nd8',...
       'consolidating q8p8 - nd24',...
       'consolidating q8p8 Alpha1 - nd8',...
       'consolidating q8p8 Alpha1 - nd24')
set(gca,'FontName','Helvetica','FontSize',16)

