
%close all
clear all

format short e

sim=load('nddataout_172.txt');
time=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

% node displ time
figure(1)
plot(time,d2,'-k','LineWidth',2)
xlabel('time (sec)')
ylabel('displacement (m)')
set(gca,'FontName','Helvetica','FontSize',16)

sim=load('nddataout_001.txt');
time=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

% node pore pressure time
figure(2)
plot(time,porepress,'-k','LineWidth',2)
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
set(gca,'FontName','Helvetica','FontSize',16)
