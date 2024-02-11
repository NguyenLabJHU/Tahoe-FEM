
%close all
clear all

format short e

sim=load('nddataout_2499_consol.txt');
time_consol=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_consol=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_2499_dyn.txt');
time_dyn=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn=sim(:,4);
d3=sim(:,5);

% node displ time
figure(1)
plot(time_consol,d2_consol,'-k',time_dyn,d2_dyn,'--k','LineWidth',2)
xlabel('time (sec)')
ylabel('displacement (m)')
legend('consolidating','dynamic')
set(gca,'FontName','Helvetica','FontSize',16)

sim=load('nddataout_2080_consol.txt');
time_consol=sim(:,1);
porepress_consol=(sim(:,2))/1000;
d1=sim(:,3);
d2_consol=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_2080_dyn.txt');
time_dyn=sim(:,1);
porepress_dyn=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn=sim(:,4);
d3=sim(:,5);

% node pore pressure time
figure(2)
plot(time_consol,porepress_consol,'-k',time_dyn,porepress_dyn,'--k','LineWidth',2)
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('consolidating','dynamic')
set(gca,'FontName','Helvetica','FontSize',16)

% figure(3)
% cosfunc=1-cos(50*time_dyn_50);
% plot(time_dyn_50,cosfunc,'-k','LineWidth',2)
% xlabel('time')
% ylabel('schedule func')
% set(gca,'FontName','Helvetica','FontSize',16)


