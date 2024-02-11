
close all
clear all

format short e

sim=load('nddata_col20m_172.txt');
time=sim(:,1);
porepress_172=(sim(:,2));
d1_172=sim(:,3);
d2_172=sim(:,4);
d3_172=sim(:,5);

sim=load('nddata_col20m_001.txt');
time=sim(:,1);
porepress_001=(sim(:,2));
d1_001=sim(:,3);
d2_001=sim(:,4);
d3_001=sim(:,5);

sim=load('nddata_col20m_352.txt');
time=sim(:,1);
porepress_352=(sim(:,2));
d1_352=sim(:,3);
d2_352=sim(:,4);
d3_352=sim(:,5);

sim=load('nddata_col20m_355.txt');
time=sim(:,1);
porepress_355=(sim(:,2));
d1_355=sim(:,3);
d2_355=sim(:,4);
d3_355=sim(:,5);

sim=load('nddata_col20m_369.txt');
time=sim(:,1);
porepress_369=(sim(:,2));
d1_369=sim(:,3);
d2_369=sim(:,4);
d3_369=sim(:,5);

% node displ time
figure(1)
plot(time,d2_001,'-k',time,d2_172,'--k',time,d2_352,'-.k','LineWidth',2)
xlabel('time (sec)')
ylabel('d2 (m)')
legend('node 1', 'node 172', 'node 352')
set(gca,'FontName','Helvetica','FontSize',16)

% node displ time
% figure(2)
% plot(time,d2_355,'-k',time,d2_369,'--k',time,d2_352,'-.k','LineWidth',2)
% xlabel('time (sec)')
% ylabel('d2 (m)')
% legend('node 355', 'node 369', 'node 352')
% set(gca,'FontName','Helvetica','FontSize',16)

% node pore pressure time
figure(3)
plot(time,porepress_001,'-k',time,porepress_172,'--k',time,porepress_352,'-.k','LineWidth',2)
xlabel('time (sec)')
ylabel('pore pressure (Pa)')
legend('node 1', 'node 172', 'node 352')
set(gca,'FontName','Helvetica','FontSize',16)

