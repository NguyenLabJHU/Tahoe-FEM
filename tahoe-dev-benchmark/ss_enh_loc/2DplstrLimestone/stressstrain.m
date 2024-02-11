clear all
%
sim=load('eldataloc_1.txt');
time=abs(sim(:,1));
sig11=abs(sim(:,2));
sig22=abs(sim(:,3));
sig12=abs(sim(:,4));
eps11=abs(sim(:,5));
eps22=abs(sim(:,6));
eps12=abs(sim(:,7));
loc_ip1=abs(sim(:,8));
loc_ip2=abs(sim(:,9));
loc_ip3=abs(sim(:,10));
%ip4
alph11=abs(sim(:,11));
alph22=abs(sim(:,12));
alph33=abs(sim(:,13));
alph23=abs(sim(:,14));
alph13=abs(sim(:,15));
alph12=abs(sim(:,16));
kappa=abs(sim(:,17));
press=abs(sim(:,18));
J2=abs(sim(:,19));
J3=abs(sim(:,20));
loc_ip4=abs(sim(:,21));
%
% figure(2)
% plot(time,alph12)
% xlabel('time')
% ylabel('backstress')
% legend('alph12')
%
% figure(3)
% plot(time,kappa)
% xlabel('time')
% ylabel('kappa')
% legend('kappa')
%
% figure(4)
% plot(time,alph11,time,alph22,time,alph33)
% xlabel('time')
% ylabel('backstress')
% legend('alph11','alph22','alph33')
%
figure(5)
plot(eps22,sig22)
xlabel('strain')
ylabel('stress')
%legend('sig22')
%
sim=load('ss_enh_isv.txt');
loc_flag=(sim(:,1));
zeta=(sim(:,2));
gamma_delta=(sim(:,3));
Q_S=(sim(:,4));
P_S=(sim(:,5));
q_St=(sim(:,6));
cohesion=(sim(:,7));
friction=(sim(:,8));
dilation=(sim(:,9));
%
figure(2)
plot(zeta,cohesion)
xlabel('zeta')
ylabel('cohesion')
%

