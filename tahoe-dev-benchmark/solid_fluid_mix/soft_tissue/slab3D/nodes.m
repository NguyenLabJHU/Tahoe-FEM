
%close all
clear all

format short e

sim=load('nddata_332.txt');
time=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3)*1000;
d2=sim(:,4)*1000;
d3=sim(:,5)*1000;

% node displ time
figure(1)
plot(time,d2,'-k','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('displacement (mm)')
legend('bottom node, center')
set(gca,'FontName','Helvetica','FontSize',16)

% node pore pressure time
figure(2)
plot(time,porepress,'-k','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('bottom node, center')
set(gca,'FontName','Helvetica','FontSize',16)
