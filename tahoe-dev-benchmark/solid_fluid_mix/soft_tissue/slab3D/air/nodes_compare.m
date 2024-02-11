
%close all
clear all

format short e

drained_time_factor = 2e-3;

%%%%center, bottom
sim=load('nddata_dyn_hold_drain_316.txt');
time_dyn_hold_drain=sim(:,1);
porepress_dyn_hold_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_hold_drain=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_dyn_drain_316.txt');
time_dyn_drain=sim(:,1);
porepress_dyn_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_drain=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_dyn_hold_316.txt');
time_dyn_hold=sim(:,1);
porepress_dyn_hold=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_hold=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_dyn_316.txt');
time_dyn=sim(:,1);
porepress_dyn=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn=sim(:,4)*1000;
d3=sim(:,5);

% node displ time
figure(1)
plot(time_dyn,d2_dyn,'-.r',time_dyn_drain,d2_dyn_drain,'-.m',time_dyn_hold,d2_dyn_hold,'-.g',time_dyn_hold_drain,d2_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('displacement (mm)')
legend('dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, center')
set(gca,'FontName','Helvetica','FontSize',16)
%axis([0 0.1 -10 2])

% node pore pressure time
figure(2)
plot(time_dyn,porepress_dyn,'-.r',time_dyn_drain,porepress_dyn_drain,'-.m',time_dyn_hold,porepress_dyn_hold,'-.g',time_dyn_hold_drain,porepress_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, center')
set(gca,'FontName','Helvetica','FontSize',16)
%axis([0 0.1 -0.3 0.4])

%%%%corner, bottom
sim=load('nddata_dyn_hold_drain_628.txt');
time_dyn_hold_drain=sim(:,1);
porepress_dyn_hold_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_hold_drain=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_dyn_drain_628.txt');
time_dyn_drain=sim(:,1);
porepress_dyn_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_drain=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_dyn_hold_628.txt');
time_dyn_hold=sim(:,1);
porepress_dyn_hold=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn_hold=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_dyn_628.txt');
time_dyn=sim(:,1);
porepress_dyn=(sim(:,2))/1000;
d1=sim(:,3);
d2_dyn=sim(:,4)*1000;
d3=sim(:,5);

% node displ time
figure(5)
plot(time_dyn,d2_dyn,'-.r',time_dyn_drain,d2_dyn_drain,'-.m',time_dyn_hold,d2_dyn_hold,'-.g',time_dyn_hold_drain,d2_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('displacement (mm)')
legend('dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, corner')
set(gca,'FontName','Helvetica','FontSize',16)
%axis([0 0.1 -10 2])

% node pore pressure time
figure(6)
plot(time_dyn,porepress_dyn,'-.r',time_dyn_drain,porepress_dyn_drain,'-.m',time_dyn_hold,porepress_dyn_hold,'-.g',time_dyn_hold_drain,porepress_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, corner')
set(gca,'FontName','Helvetica','FontSize',16)
%axis([0 0.1 -0.3 0.4])

