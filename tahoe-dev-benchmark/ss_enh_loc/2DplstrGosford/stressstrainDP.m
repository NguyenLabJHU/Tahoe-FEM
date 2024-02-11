clear all
%
sim=load('eldata_1.txt');
time=(sim(:,1));
%ip1
sig11_1=(sim(:,2));
sig22_1=(sim(:,3));
sig12_1=(sim(:,4));
eps11_1=(sim(:,5));
eps22_1=(sim(:,6));
eps12_1=(sim(:,7));
%ip1
alpha_1=(sim(:,8));
plstr_1=(sim(:,9));
VM_1=(sim(:,10));
press_1=(sim(:,11));
%
% initial yield surface
alpha0=1.3e4;
alpha=alpha0;
beta=0.5; %radians, =30degrees
p=-1e5:1:2.5e4;
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
