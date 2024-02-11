clear all
%
sim=load('eldataMRenh_1.txt');
time=(sim(:,1));
%ip1
sig11_1=(sim(:,2));
sig22_1=(sim(:,3));
sig12_1=(sim(:,4));
eps11_1=(sim(:,5));
eps22_1=(sim(:,6));
eps12_1=(sim(:,7));
%ip1
tension_1=(sim(:,8));
cohesion_1=(sim(:,9));
friction_1=(sim(:,10));
dilation_1=(sim(:,11));
VM_1=(sim(:,12));
press_1=(sim(:,13));
%
% initial yield surface
c=1.3e4;
chi=2.25e4;
phi=0.5236; %radians, =30degrees
p=-1e5:1:2.5e4;
tau2=(c-p*tan(phi)).^2 - (c-chi*tan(phi))^2;
tau2pos=0.5*(tau2+abs(tau2));
taupos0=sqrt(tau2pos);
%
% final yield surface
c=cohesion_1(length(cohesion_1));
chi=tension_1(length(cohesion_1));
tanphi=friction_1(length(cohesion_1)); 
tau2=(c-p*tanphi).^2 - (c-chi*tanphi)^2;
tau2pos=0.5*(tau2+abs(tau2));
taupos1=sqrt(tau2pos);
%
factor=1/sqrt(3);
%
figure(1)
plot(press_1,VM_1*factor,'o',p,taupos0,p,taupos1,'LineWidth',1)
xlabel('mean stress (MPa)')
ylabel('sqrt(J2) (MPa)')
legend('stress path','initial MR yield surface','final MR yield surface')
set(gca,'FontName','Helvetica','FontSize',16)
%
