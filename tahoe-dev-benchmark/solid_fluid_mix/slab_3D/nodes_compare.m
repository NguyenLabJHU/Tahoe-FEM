
%close all
clear all

format short e

sim=load('nddata_331_consol.txt');
time_consol=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_consol=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_331_drain.txt');
time_drain=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2_drain=sim(:,4)*1000;
d3=sim(:,5);

% node displ time
figure(1)
plot(time_drain,d2_drain,'-k',time_consol,d2_consol,'--k','LineWidth',2)
xlabel('time (sec)')
ylabel('displacement (mm)')
legend('drained','consolidating')
set(gca,'FontName','Helvetica','FontSize',16)

sim=load('nddata_328_consol.txt');
time_consol=sim(:,1);
porepress_consol=(sim(:,2))/1000;
d1=sim(:,3);
d2_consol=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_328_drain.txt');
time_drain=sim(:,1);
porepress_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_drain=sim(:,4)*1000;
d3=sim(:,5);

% node pore pressure time
figure(2)
plot(time_drain,porepress_drain,'-k',time_consol,porepress_consol,'--k','LineWidth',2)
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('drained','consolidating')
set(gca,'FontName','Helvetica','FontSize',16)

