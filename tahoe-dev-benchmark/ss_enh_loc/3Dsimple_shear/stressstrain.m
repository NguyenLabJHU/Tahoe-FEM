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
plstr=(sim(:,15));
VM=(sim(:,16));
press=(sim(:,17));
loccheck=(sim(:,18));
ellocflag=(sim(:,19));
%
% figure(1)
% %plot(abs(eps12),VM)
% plot(abs(eps12),sig12)
% xlabel('strain')
% ylabel('shear stress')
%ylabel('VM stress')
%legend('100 steps','40 steps')
%
sim=load('nddata_7.txt');
time=(sim(:,1));
displ=(sim(:,2));
%equal force for 4 nodes
force=4*(sim(:,3));
%
figure(2)
plot(displ,force)
xlabel('displ')
ylabel('force')
%
sim=load('ss_enh_isv.txt');
elem=(sim(:,1));
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
figure(3)
plot(zeta,cohesion)
xlabel('zeta')
ylabel('cohesion')
%
