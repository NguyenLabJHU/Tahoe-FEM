clear all
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
figure(8)
plot(d,f)
xlabel('DISPLACEMENT (mm)') 
ylabel('FORCE (MN)')
%legend('node 27')
%
sim=load('ss_enh_isv.txt');
loc_flag=(sim(:,1));
zeta=(sim(:,2));
gamma_delta=(sim(:,3));
Q_S=(sim(:,4));
P_S=(sim(:,5));
q_St=(sim(:,6));
cohesion=(sim(:,7));
friction=(sim(:,8));
dilation=(sim(:,9));
%
figure(2)
plot(zeta,cohesion)
xlabel('ZETA (m)')
ylabel('COHESION (MPa)')
%
figure(3)
plot(zeta,friction)
xlabel('ZETA (m)')
ylabel('FRICTION ANGLE (radian)')
%
figure(4)
plot(zeta,dilation)
xlabel('ZETA (m)')
ylabel('DILATION ANGLE (radian)')
%


