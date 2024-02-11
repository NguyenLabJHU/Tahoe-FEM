clear all
%

%
sim=load('nddata_27.txt');
time=abs(sim(:,1));
dx=abs(sim(:,2));
dy=abs(sim(:,3));
dz=abs(sim(:,4));
d=1000*sqrt(dx.^2 + dy.^2 + dz.^2);
fx=abs(sim(:,5));
fy=abs(sim(:,6));
fz=abs(sim(:,7));
f=sqrt(fx.^2 + fy.^2 + fz.^2);
%
figure(7)
plot(time,f,'LineWidth',1)
xlabel('TIME (sec)')
%set(get(gco,'XLabel'),'DISPLACEMENT (mm)','FontName','Helvetica','FontSize',16)
ylabel('FORCE (MN)')
%legend('node 27')
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(8)
plot(d,f,'LineWidth',1)
xlabel('DISPLACEMENT (mm)')
%set(get(gco,'XLabel'),'DISPLACEMENT (mm)','FontName','Helvetica','FontSize',16)
ylabel('FORCE (MN)')
%legend('node 27')
set(gca,'FontName','Helvetica','FontSize',16)
%
sim=load('ss_enh_isv.txt');
elem=(sim(:,1));
loc_flag=(sim(:,2));
zeta=(sim(:,3));
gamma_delta=(sim(:,4));
Q_S=(sim(:,5));
P_S=(sim(:,6));
q_St=(sim(:,7));
cohesion=(sim(:,8));
friction=(sim(:,9));
dilation=(sim(:,10));
%
figure(2)
plot(zeta,cohesion,'LineWidth',1)
xlabel('ZETA (m)')
ylabel('COHESION (MPa)')
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(3)
plot(zeta,friction,'LineWidth',1)
xlabel('ZETA (m)')
ylabel('FRICTION ANGLE (radian)')
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(4)
plot(zeta,dilation,'LineWidth',1)
xlabel('ZETA (m)')
ylabel('DILATION ANGLE (radian)')
set(gca,'FontName','Helvetica','FontSize',16)
%



