clear all
%
sim=load('eldataMRenh_1.txt');
time=abs(sim(:,1));
sig11=abs(sim(:,2));
sig22=abs(sim(:,3));
sig12=abs(sim(:,4));
eps11=abs(sim(:,5));
eps22=abs(sim(:,6));
eps12=abs(sim(:,7));
chi=abs(sim(:,8));
cohesion=abs(sim(:,9));
friction=abs(sim(:,10));
dilation=abs(sim(:,11));
VM=abs(sim(:,12));
press=abs(sim(:,13));
%
figure(10)
plot(time,sig22)
xlabel('time')
ylabel('stress')
%legend('sig22')
%
% figure(5)
% plot(eps22,sig22)
% xlabel('strain')
% ylabel('stress')
% %legend('sig22')
%
sim=load('nddataMRenh_4.txt');
time=(sim(:,1));
%node 4
dy=(sim(:,2));
fy=(sim(:,3));
%
figure(11)
plot(-dy,-2*fy)
xlabel('displ')
ylabel('force')
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
figure(12)
plot(zeta,cohesion)
xlabel('zeta (m)')
ylabel('cohesion (Pa)')
%
figure(13)
plot(zeta,friction)
xlabel('zeta (m)')
ylabel('friction (rad)')
%
figure(14)
plot(zeta,dilation)
xlabel('zeta (m)')
ylabel('dilation (rad) ')
%

