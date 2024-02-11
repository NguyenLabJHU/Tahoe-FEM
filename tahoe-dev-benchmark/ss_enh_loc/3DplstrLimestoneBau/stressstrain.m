clear all
%
sim=load('eldataloc_1_bau.txt');
sig22_bau=(sim(:,3));
eps22_bau=(sim(:,9));
%
sim=load('eldataloc_1_mono.txt');
sig22_mono=(sim(:,3));
eps22_mono=(sim(:,9));
%
figure(1)
plot(-eps22_mono,-sig22_mono,-eps22_bau,-sig22_bau,'LineWidth',2)
xlabel('STRAIN')
ylabel('STRESS (MPa)')
legend('monotonic loading','reverse loading')
set(gca,'FontName','Helvetica','FontSize',16)
%
sim=load('eldataloc_1_bau_j3.txt');
sig22_bau_j3=(sim(:,3));
eps22_bau_j3=(sim(:,9));
%
sim=load('eldataloc_1_mono_j3.txt');
sig22_mono_j3=(sim(:,3));
eps22_mono_j3=(sim(:,9));
%
figure(2)
plot(-eps22_mono,-sig22_mono,-eps22_mono_j3,-sig22_mono_j3,-eps22_bau_j3,-sig22_bau_j3,'LineWidth',2)
xlabel('STRAIN')
ylabel('STRESS (MPa)')
legend('monotonic loading','monotonic loading - J3','reverse loading - J3')
set(gca,'FontName','Helvetica','FontSize',16)
%