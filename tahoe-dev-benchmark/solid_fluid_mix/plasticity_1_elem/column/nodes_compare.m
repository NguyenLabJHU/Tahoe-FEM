
%close all
clear all

format short e

sim=load('nddataout_172_consol_elastic.txt');
time=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_elastic=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_172_consol.txt');
time=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_172_dyn_pi.txt');
time_dyn_pi=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_pi=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_172_dyn_50.txt');
time_dyn_50=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_50=sim(:,4);
d3=sim(:,5);

% node displ time
figure(1)
plot(time,d2_elastic,'-k',time,d2,'--k',time_dyn_pi,d2_dyn_pi,'-.k',time_dyn_50,d2_dyn_50,':k','LineWidth',2)
xlabel('time')
ylabel('d2')
legend('elastic','plastic','plastic dynamic Pi','plastic dynamic 50')
set(gca,'FontName','Helvetica','FontSize',16)

sim=load('nddataout_1_consol_elastic.txt');
time=sim(:,1);
porepress_elastic=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_1_consol.txt');
time=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_1_dyn_pi.txt');
time_dyn_pi=sim(:,1);
porepress_dyn_pi=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_1_dyn_50.txt');
time_dyn_50=sim(:,1);
porepress_dyn_50=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

% node pore pressure time
figure(2)
plot(time,porepress_elastic,'-k',time,porepress,'--k',time_dyn_pi,porepress_dyn_pi,'-.k',time_dyn_50,porepress_dyn_50,':k','LineWidth',2)
xlabel('time')
ylabel('pore pressure')
legend('elastic','plastic','plastic dynamic Pi','plastic dynamic 50')
set(gca,'FontName','Helvetica','FontSize',16)



% figure(3)
% cosfunc=1-cos(50*time_dyn_50);
% plot(time_dyn_50,cosfunc,'-k','LineWidth',2)
% xlabel('time')
% ylabel('schedule func')
% set(gca,'FontName','Helvetica','FontSize',16)


