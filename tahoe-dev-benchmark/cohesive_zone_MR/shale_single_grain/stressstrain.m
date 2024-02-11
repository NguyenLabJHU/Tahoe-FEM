clear all
%
sim=load('eldata_812.txt');
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
tension_1=(sim(:,26));
cohesion_1=(sim(:,27));
friction_1=(sim(:,28));
dilation_1=(sim(:,29));
VM_1=(sim(:,30));
press_1=(sim(:,31));
%ip2
tension_2=(sim(:,32));
cohesion_2=(sim(:,33));
friction_2=(sim(:,34));
dilation_2=(sim(:,35));
VM_2=(sim(:,36));
press_2=(sim(:,37));
%ip3
tension_3=(sim(:,38));
cohesion_3=(sim(:,39));
friction_3=(sim(:,40));
dilation_3=(sim(:,41));
VM_3=(sim(:,42));
press_3=(sim(:,43));
%ip4
tension_4=(sim(:,44));
cohesion_4=(sim(:,45));
friction_4=(sim(:,46));
dilation_4=(sim(:,47));
VM_4=(sim(:,48));
press_4=(sim(:,49));
%
% initial yield surface
c=5;
chi=1;
phi=0.436; %radians, =25degrees
p=-15:.01:1.1;
tau2=(c-p*tan(phi)).^2 - (c-chi*tan(phi))^2;
tau2pos=0.5*(tau2+abs(tau2));
taupos0=sqrt(tau2pos);
%
% final yield surface
c=cohesion_4(length(cohesion_4));
chi=tension_4(length(cohesion_4));
tanphi=friction_4(length(cohesion_4)); 
p=-15:.01:1.1;
tau2=(c-p*tanphi).^2 - (c-chi*tanphi)^2;
tau2pos=0.5*(tau2+abs(tau2));
taupos1=sqrt(tau2pos);
%
factor=1/sqrt(3);
%
figure(1)
plot(press_4,VM_4*factor,'ok',p,taupos0,'--k',p,taupos1,'-k','LineWidth',1)
xlabel('mean stress (MPa)')
ylabel('sqrt(J2) (MPa)')
legend('stress path','initial yield surface','final yield surface')
set(gca,'FontName','Helvetica','FontSize',16)
%
