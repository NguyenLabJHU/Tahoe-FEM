clear all
%
sim=load('eldata_25.txt');
time=(sim(:,1));
%ip1
sig11_1=(sim(:,2));
sig22_1=(sim(:,3));
sig12_1=(sim(:,4));
eps11_1=(sim(:,5));
eps22_1=(sim(:,6));
eps12_1=(sim(:,7));
%ip2
sig11_2=(sim(:,8));
sig22_2=(sim(:,9));
sig12_2=(sim(:,10));
eps11_2=(sim(:,11));
eps22_2=(sim(:,12));
eps12_2=(sim(:,13));
%ip3
sig11_3=(sim(:,14));
sig22_3=(sim(:,15));
sig12_3=(sim(:,16));
eps11_3=(sim(:,17));
eps22_3=(sim(:,18));
eps12_3=(sim(:,19));
%ip4
sig11_4=(sim(:,20));
sig22_4=(sim(:,21));
sig12_4=(sim(:,22));
eps11_4=(sim(:,23));
eps22_4=(sim(:,24));
eps12_4=(sim(:,25));
%ip1
alpha_1=(sim(:,26));
VM_1=(sim(:,27));
press_1=(sim(:,28));
%
% initial yield surface
alpha0=6.1;
alpha=alpha0;
beta=0.436; %radians, =30degrees
p=-1e1:1:2.5e1;
tau=sqrt(3)*(alpha-beta*p);
taupos0=0.5*(tau+abs(tau));
%
% final yield surface
alpha1=alpha_1(length(alpha_1))/sqrt(3);
alpha=alpha0-alpha1;
tau=sqrt(3)*(alpha-beta*p);
taupos1=0.5*(tau+abs(tau));
%
factor=1/sqrt(3);
%
figure(2)
plot(press_1,VM_1*factor,'o',p,taupos0*factor,p,taupos1*factor,'LineWidth',1)
%plot(press_1,VM_1,'o',p,taupos0,p,taupos1,'LineWidth',1)
xlabel('mean stress (MPa)')
ylabel('sqrt(J2) (MPa)')
%ylabel('VM (MPa)')
legend('stress path','initial DP yield surface','final DP yield surface')
set(gca,'FontName','Helvetica','FontSize',16)
%
