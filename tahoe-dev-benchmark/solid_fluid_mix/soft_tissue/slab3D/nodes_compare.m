
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

sim=load('nddata_consol_316.txt');
time_consol=sim(:,1);
porepress_consol=(sim(:,2))/1000;
d1=sim(:,3);
d2_consol=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_drain_316.txt');
time_drain=sim(:,1)*drained_time_factor;
porepress_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_drain=sim(:,4)*1000;
d3=sim(:,5);

% node displ time
figure(1)
plot(time_drain,d2_drain,'-k',time_consol,d2_consol,'--b',time_dyn,d2_dyn,'-.r',time_dyn_drain,d2_dyn_drain,'-.m',time_dyn_hold,d2_dyn_hold,'-.g',time_dyn_hold_drain,d2_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('displacement (mm)')
legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, center')
set(gca,'FontName','Helvetica','FontSize',16)
axis([0 0.1 -10 2])

% node pore pressure time
figure(2)
plot(time_drain,porepress_drain,'-k',time_consol,porepress_consol,'--b',time_dyn,porepress_dyn,'-.r',time_dyn_drain,porepress_dyn_drain,'-.m',time_dyn_hold,porepress_dyn_hold,'-.g',time_dyn_hold_drain,porepress_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, center')
set(gca,'FontName','Helvetica','FontSize',16)
axis([0 0.1 -0.3 0.4])

%%%%center, top
% sim=load('nddata_dyn_drain_317.txt');
% time_dyn_drain=sim(:,1);
% porepress_dyn_drain=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_dyn_drain=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_dyn_hold_317.txt');
% time_dyn_hold=sim(:,1);
% porepress_dyn_hold=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_dyn_hold=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_dyn_317.txt');
% time_dyn=sim(:,1);
% porepress_dyn=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_dyn=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_consol_317.txt');
% time_consol=sim(:,1);
% porepress_consol=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_consol=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_drain_317.txt');
% time_drain=sim(:,1)*drained_time_factor;
% porepress_drain=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_drain=sim(:,4)*1000;
% d3=sim(:,5);
% 
% % node displ time
% figure(3)
% plot(time_drain,d2_drain,'-k',time_consol,d2_consol,'--b',time_dyn,d2_dyn,'-.r',time_dyn_drain,d2_dyn_drain,'-.m',time_dyn_hold,d2_dyn_hold,'-.g',time_dyn_hold_drain,d2_dyn_hold_drain,'-.c','LineWidth',2)
% grid on
% xlabel('time (sec)')
% ylabel('displacement (mm)')
% legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
% title('top node, center')
% set(gca,'FontName','Helvetica','FontSize',16)
% 
% % node pore pressure time
% figure(4)
% plot(time_drain,porepress_drain,'-k',time_consol,porepress_consol,'--b',time_dyn,porepress_dyn,'-.r',time_dyn_drain,porepress_dyn_drain,'-.m',time_dyn_hold,porepress_dyn_hold,'-.g',time_dyn_hold_drain,porepress_dyn_hold_drain,'-.c','LineWidth',2)
% grid on
% xlabel('time (sec)')
% ylabel('pore pressure (kPa)')
% legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
% title('top node, center')
% set(gca,'FontName','Helvetica','FontSize',16)

%%%%corner, bottom
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

sim=load('nddata_consol_628.txt');
time_consol=sim(:,1);
porepress_consol=(sim(:,2))/1000;
d1=sim(:,3);
d2_consol=sim(:,4)*1000;
d3=sim(:,5);

sim=load('nddata_drain_628.txt');
time_drain=sim(:,1)*drained_time_factor;
porepress_drain=(sim(:,2))/1000;
d1=sim(:,3);
d2_drain=sim(:,4)*1000;
d3=sim(:,5);

% node displ time
figure(5)
plot(time_drain,d2_drain,'-k',time_consol,d2_consol,'--b',time_dyn,d2_dyn,'-.r',time_dyn_drain,d2_dyn_drain,'-.m',time_dyn_hold,d2_dyn_hold,'-.g',time_dyn_hold_drain,d2_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('displacement (mm)')
legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, corner')
set(gca,'FontName','Helvetica','FontSize',16)
axis([0 0.1 -10 2])

% node pore pressure time
figure(6)
plot(time_drain,porepress_drain,'-k',time_consol,porepress_consol,'--b',time_dyn,porepress_dyn,'-.r',time_dyn_drain,porepress_dyn_drain,'-.m',time_dyn_hold,porepress_dyn_hold,'-.g',time_dyn_hold_drain,porepress_dyn_hold_drain,'-.c','LineWidth',2)
grid on
xlabel('time (sec)')
ylabel('pore pressure (kPa)')
legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
title('bottom node, corner')
set(gca,'FontName','Helvetica','FontSize',16)
axis([0 0.1 -0.3 0.4])

%%%%corner, top
% sim=load('nddata_dyn_drain_631.txt');
% time_dyn_drain=sim(:,1);
% porepress_dyn_drain=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_dyn_drain=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_dyn_hold_631.txt');
% time_dyn_hold=sim(:,1);
% porepress_dyn_hold=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_dyn_hold=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_dyn_631.txt');
% time_dyn=sim(:,1);
% porepress_dyn=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_dyn=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_consol_631.txt');
% time_consol=sim(:,1);
% porepress_consol=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_consol=sim(:,4)*1000;
% d3=sim(:,5);
% 
% sim=load('nddata_drain_631.txt');
% time_drain=sim(:,1)*drained_time_factor;
% porepress_drain=(sim(:,2))/1000;
% d1=sim(:,3);
% d2_drain=sim(:,4)*1000;
% d3=sim(:,5);
% 
% % node displ time
% figure(7)
% plot(time_drain,d2_drain,'-k',time_consol,d2_consol,'--b',time_dyn,d2_dyn,'-.r',time_dyn_drain,d2_dyn_drain,'-.m',time_dyn_hold,d2_dyn_hold,'-.g',time_dyn_hold_drain,d2_dyn_hold_drain,'-.c','LineWidth',2)
% grid on
% xlabel('time (sec)')
% ylabel('displacement (mm)')
% legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
% title('top node, corner')
% set(gca,'FontName','Helvetica','FontSize',16)
% 
% % node pore pressure time
% figure(8)
% plot(time_drain,porepress_drain,'-k',time_consol,porepress_consol,'--b',time_dyn,porepress_dyn,'-.r',time_dyn_drain,porepress_dyn_drain,'-.m',time_dyn_hold,porepress_dyn_hold,'-.g',time_dyn_hold_drain,porepress_dyn_hold_drain,'-.c','LineWidth',2)
% grid on
% xlabel('time (sec)')
% ylabel('pore pressure (kPa)')
% legend('drained','consolidating','dynamic impulse','dynamic impulse drained','dynamic hold','dynamic hold drained')
% title('top node, corner')
% set(gca,'FontName','Helvetica','FontSize',16)

