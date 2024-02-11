clear all
%
sim=load('eldata_1.txt');
time=(sim(:,1));
%ip1
sig11=(sim(:,2));
sig22=(sim(:,3));
sig33=(sim(:,4));
sig23=(sim(:,5));
sig13=(sim(:,6));
sig12=(sim(:,7));
eps11=(sim(:,8));
eps22=(sim(:,9));
eps33=(sim(:,10));
eps23=(sim(:,11));
eps13=(sim(:,12));
eps12=(sim(:,13));
kappa=(sim(:,14));
VM=(sim(:,15));
press=(sim(:,16));
loccheck=(sim(:,17));
%
% figure(1)
% plot(abs(eps22),abs(sig22))
% xlabel('strain')
% ylabel('stress')
%
sim=load('ss_enh_isv.txt');
elem_num=(sim(:,1));
loc_flag=(sim(:,2));
zeta=(sim(:,3));
gamma_delta=(sim(:,4));
Q_S=(sim(:,5));
P_S=(sim(:,6));
q_St=(sim(:,7));
cohesion=(sim(:,8));
friction=(sim(:,9));
dilation=(sim(:,10));
%
sim=load('./cse_comparison/ep/nddata_isv_3_1e9.txt');
time_cse_ep=(sim(:,1));
%node 3
upt_cse_ep=(sim(:,2));
upn_cse_ep=(sim(:,3));
up_mag_cse_ep9=sqrt(upt_cse_ep.^2+upn_cse_ep.^2);
chi_cse_ep=(sim(:,4));
cohesion_cse_ep9=(sim(:,5));
friction_cse_ep9=(sim(:,6));
dilation_cse_ep9=(sim(:,7));
%
sim=load('./cse_comparison/ep/nddata_isv_3_1e8.txt');
time_cse_ep=(sim(:,1));
%node 3
upt_cse_ep=(sim(:,2));
upn_cse_ep=(sim(:,3));
up_mag_cse_ep8=sqrt(upt_cse_ep.^2+upn_cse_ep.^2);
chi_cse_ep=(sim(:,4));
cohesion_cse_ep8=(sim(:,5));
friction_cse_ep8=(sim(:,6));
dilation_cse_ep8=(sim(:,7));
%
up_mag_cse_rp=load('./rp_matlab_cse/rp_cse_upmag.txt');
chi_cse_rp=load('./rp_matlab_cse/rp_cse_chi.txt');
cohesion_cse_rp=load('./rp_matlab_cse/rp_cse_c.txt');
friction_cse_rp=load('./rp_matlab_cse/rp_cse_phi.txt');
dilation_cse_rp=load('./rp_matlab_cse/rp_cse_psi.txt');
%
figure(2)
plot(zeta,cohesion,'-k',up_mag_cse_ep9,cohesion_cse_ep9,'--r',up_mag_cse_ep8,cohesion_cse_ep8,'-.r',...
    up_mag_cse_rp,cohesion_cse_rp,'ob','LineWidth',1.5)
xlabel('PLASTIC JUMP DISPLACEMENT MAGNITUDE (mm)')
ylabel('COHESION (Pa)')
legend('embedded discontinuity','cse EP: k_n = k_t = 1e9 N/m','cse EP: k_n = k_t = 1e8 N/m','cse RP')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','0.5','1','1.5','2','2.5','3'})
set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(3)
plot(zeta,friction,'-k',up_mag_cse_ep9,friction_cse_ep9,'--r',up_mag_cse_ep8,friction_cse_ep8,'-.r',...
    up_mag_cse_rp,friction_cse_rp,'ob','LineWidth',1.5)
xlabel('PLASTIC JUMP DISPLACEMENT MAGNITUDE (mm)')
ylabel('FRICTION ANGLE (rad)')
legend('embedded discontinuity','cse EP: k_n = k_t = 1e9 N/m','cse EP: k_n = k_t = 1e8 N/m','cse RP')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','0.5','1','1.5','2','2.5','3'})
set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(4)
plot(zeta,dilation,'-k',up_mag_cse_ep9,dilation_cse_ep9,'--r',up_mag_cse_ep8,dilation_cse_ep8,'-.r',...
    up_mag_cse_rp,dilation_cse_rp,'ob','LineWidth',1.5)
xlabel('PLASTIC JUMP DISPLACEMENT MAGNITUDE (mm)')
ylabel('DILATION ANGLE (rad)')
legend('embedded discontinuity','cse EP: k_n = k_t = 1e9 N/m','cse EP: k_n = k_t = 1e8 N/m','cse RP')
%axis([0 2.5 0 30])
set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
%
sim=load('nddata_2.txt');
time=(sim(:,1));
%node 2
dy=(sim(:,2));
fy=(sim(:,3));
%
sim=load('./cse_comparison/ep/nddata_7_1e9.txt');
time_cse_ep9=(sim(:,1));
%node 7
dy_cse_ep9=(sim(:,2));
fy_cse_ep9=(sim(:,3));
%
sim=load('./cse_comparison/ep/nddata_7_1e8.txt');
time_cse_ep8=(sim(:,1));
%node 7
dy_cse_ep8=(sim(:,2));
fy_cse_ep8=(sim(:,3));
%
d_rp=load('./rp_matlab_cse/rp_cse_d.txt');
F_rp=load('./rp_matlab_cse/rp_cse_F.txt');
%
figure(5)
%there are four nodes, so multiply by 4 for total vertical reaction force
plot(abs(dy),4*abs(fy)/0.08,'-k',abs(dy_cse_ep9),2*abs(fy_cse_ep9),'--r',...
    abs(dy_cse_ep8),2*abs(fy_cse_ep8),'-.r',d_rp,F_rp,'ob','LineWidth',1.5)
xlabel('AXIAL DISPLACEMENT (mm)')
ylabel('AXIAL FORCE (N/m)')
legend('embedded discontinuity','cse EP: k_n = k_t = 1e9 N/m','cse EP: k_n = k_t = 1e8 N/m','cse RP')
%axis([0 2.5 0 30])
set(gca,'XTickLabel',{'0','1','2','3','4','5'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
%

