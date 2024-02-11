clear all
%
data=load('ta18.dat');
stressexp=(1e6)*data(:,1);
strainexp=data(:,2);
%
% sim=load('nddata_vel_27.txt');
% time=(sim(:,1))/1000;
% d1=abs(sim(:,2));
% d2=abs(sim(:,3));
% d3=abs(sim(:,4));
% sig11=abs(sim(:,5));
% sig22=abs(sim(:,6));
% sig33=abs(sim(:,7));
% sig23=(sim(:,8));
% sig13=(sim(:,9));
% sig12=(sim(:,10));
% kappa=(sim(:,11));
% phi=(sim(:,12));
% eqps=(sim(:,13));
% temp=(sim(:,14));
% %
% figure(1)
% plot(strainexp,stressexp,'x',time,sig33,time,kappa)
% xlabel('strain')
% ylabel('stress')
% %legend('data','sig11','sig22','sig33','kappa')
% legend('data','sig33','kappa')
% %
% figure(2)
% plot(time,temp)
% xlabel('time')
% ylabel('temp')
% legend('temp')
% %
sim=load('eldata_1.txt');
time=(sim(:,1))/1000;
sig11=abs(sim(:,2));
sig22=abs(sim(:,3));
sig33=abs(sim(:,4));
sig23=(sim(:,5));
sig13=(sim(:,6));
sig12=(sim(:,7));
eps11=abs(sim(:,8));
eps22=abs(sim(:,9));
eps33=abs(sim(:,10));
eps23=(sim(:,11));
eps13=(sim(:,12));
eps12=(sim(:,13));
kappa=(sim(:,14));
phi=(sim(:,15));
eqps=(sim(:,16));
temp=(sim(:,17));
%
% sim=load('eldata_1_old.txt');
% time_old=(sim(:,1))/1000;
% sig33_old=abs(sim(:,4));
%
figure(3)
%plot(strainexp,stressexp,'x',eps11,sig11,eps22,sig22,eps33,sig33,eps33,kappa)
plot(strainexp,stressexp,'x',time,sig33,time,kappa)
%plot(strainexp,stressexp,'x',time,sig33,time_old,sig33_old,time,kappa)
xlabel('strain')
ylabel('stress')
legend('data','sig33','kappa')
%legend('data','sig33', 'sig33-old', 'kappa')
%
figure(4)
plot(time,temp)
xlabel('strain')
ylabel('temp')
legend('temp')

