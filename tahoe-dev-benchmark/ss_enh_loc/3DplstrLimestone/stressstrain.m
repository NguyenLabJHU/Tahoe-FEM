clear all
%
sim=load('eldataloc_1_2pt5.txt');
sig22_2pt5=(sim(:,3));
eps22_2pt5=(sim(:,9));
%
sim=load('eldataloc_1_assoc.txt');
sig22_assoc=(sim(:,3));
eps22_assoc=(sim(:,9));
%
%L=1e-4 1/MPa
sim=load('eldataloc_1_pt25.txt');
sig22_pt25=(sim(:,3));
eps22_pt25=(sim(:,9));
%
sim=load('eldataloc_1_pt025.txt');
sig22_pt025=(sim(:,3));
eps22_pt025=(sim(:,9));
%
sim=load('eldataloc_1_nonassoc_weak.txt');
sig22_nonassoc=(sim(:,3));
eps22_nonassoc=(sim(:,9));
%
sim=load('eldataloc_1_nonassoc_strong.txt');
sig22_nonassoc_st=(sim(:,3));
eps22_nonassoc_st=(sim(:,9));
%
figure(1)
plot(abs(eps22_pt025),abs(sig22_pt025),abs(eps22_pt25),abs(sig22_pt25),abs(eps22_2pt5),abs(sig22_2pt5))
xlabel('strain')
ylabel('stress')
legend('0.025/sec','0.25/sec','2.5/sec')
%
figure(2)
plot(abs(eps22_assoc),abs(sig22_assoc),abs(eps22_pt25),abs(sig22_pt25),abs(eps22_pt025),abs(sig22_pt025),abs(eps22_nonassoc),abs(sig22_nonassoc),'LineWidth',2)
xlabel('STRAIN')
ylabel('STRESS (MPa)')
legend('inviscid - associative','0.25/sec - nonassociative','0.025/sec - nonassociative','inviscid - nonassociative')
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(3)
plot(abs(eps22_nonassoc),abs(sig22_nonassoc),abs(eps22_nonassoc_st),abs(sig22_nonassoc_st),'LineWidth',2)
xlabel('STRAIN')
ylabel('STRESS (MPa)')
legend('inviscid - nonassociative - weak','inviscid - nonassociative - strong')
set(gca,'FontName','Helvetica','FontSize',16)
%