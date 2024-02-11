clear all
%
sim=load('eldataloc_1.txt');
%end_step=38;
end_step=length(sim(:,1))
time=(sim(1:end_step,1));
sig11=(sim(1:end_step,2));
sig22=(sim(1:end_step,3));
sig33=(sim(1:end_step,4));
sig23=(sim(1:end_step,5));
sig13=(sim(1:end_step,6));
sig12=(sim(1:end_step,7));
eps11=(sim(1:end_step,8));
eps22=(sim(1:end_step,9));
eps33=(sim(1:end_step,10));
eps23=(sim(1:end_step,11));
eps13=(sim(1:end_step,12));
eps12=(sim(1:end_step,13));
%ip1
alph11=(sim(1:end_step,14));
alph22=(sim(1:end_step,15));
alph33=(sim(1:end_step,16));
alph23=(sim(1:end_step,17));
alph13=(sim(1:end_step,18));
alph12=(sim(1:end_step,19));
kappa=(sim(1:end_step,20));
press=sim(1:end_step,21);
J2=(sim(1:end_step,22));
J3=(sim(1:end_step,23));
locflag=(sim(:,24));
%
I1s=3*press;
devsig11=sig11-press;
devsig22=sig22-press;
devsig33=sig33-press;
xi11=devsig11-alph11;
xi22=devsig22-alph22;
xi33=devsig33-alph33;
xi23=sig23-alph23;
xi12=sig12-alph12;
xi13=sig13-alph13;
xiprod=xi11.^2+xi22.^2+xi33.^2+2*xi12.^2+2*xi23.^2+2*xi13.^2;
J2xi=.5*xiprod;
sJ2xi=sqrt(J2xi);
sJ2=sqrt(J2);
%
A=689.2;  % MPa
B=3.94e-4; % 1/MPa
C=675.2;  %MPa
R=28;    % unitless
theta=0; % radians
N=6.0;     % MPa
D1=1.47e-3;
D2=0.0;
W=.08;

%X0=-460; % MPa
kappa0=-8.05; %MPa
I1=-500:1:50;
%shear failure surface
Ff=A-C*exp(B*I1)-theta*I1;
%Ff=A-C*exp(B*I1);
Ff=0.5*(abs(Ff)+Ff);
%cap surface
X0=kappa0-R*(A-C*exp(B*kappa0)-theta*kappa0)
Fc0=1+(I1-kappa0).*(abs(I1-kappa0)-(I1-kappa0))/(2*(X0-kappa0)^2);
Fc0=0.5*(abs(Fc0)+Fc0);
kappa_last=(sim(length(time),20));
X=kappa_last-R*(A-C*exp(B*kappa_last)-theta*kappa_last)
Fc=1+(I1-kappa_last).*(abs(I1-kappa_last)-(I1-kappa_last))/(2*(X-kappa_last)^2);
Fc=0.5*(abs(Fc)+Fc);
Fy=0.5*(abs(Ff-N)+Ff-N);
fJ2yield0=Fc0.*Fy.^2;
fJ2yieldsq0=sqrt(Fc0).*Fy;
fJ2failsq0=sqrt(Fc0).*Ff;
fJ2yieldsq=sqrt(Fc).*Fy;
fJ2failsq=sqrt(Fc).*Ff;
%
% figure(2)
% plot(time,alph12)
% xlabel('strain')
% ylabel('backstress')
% legend('alph12')
% %
% figure(3)
% plot(time,kappa)
% xlabel('time')
% ylabel('kappa')
% legend('kappa')
% %
% figure(4)
% plot(time,alph11,time,alph22,time,alph33)
% xlabel('time')
% ylabel('backstress')
% legend('alph11','alph22','alph33')
% %
figure(5)
plot(-eps22,-sig22,'LineWidth',1)
xlabel('STRAIN')
ylabel('STRESS (MPa)')
%legend('100 steps','40 steps')
set(gca,'FontName','Helvetica','FontSize',16)
%
% figure(6)
% plot(time,sig22)
% xlabel('time')
% ylabel('stress')
% legend('sig22')
%
% figure(7)
% plot(time,loc_ip8,'o')
% xlabel('time')
% ylabel('locflag')
% legend('loc-ip8')
% axis([0 1 -.1 1.1])
%
% figure(8)
% plot(I1s,sJ2,'o',I1s,sJ2xi)
% xlabel('I1')
% ylabel('sJ2')
% legend('sJ2','sJ2xi')
%
% figure(8)
% plot(I1s,sJ2,'o',I1,fJ2yieldsq0,I1,fJ2failsq0,'--','LineWidth',1)
% xlabel('I1')
% ylabel('sJ2')
% legend('stress path','initial yield surface','initial failure surface')
% set(gca,'FontName','Helvetica','FontSize',16)
%
figure(9)
%plot(I1s,sJ2xi,'o',I1,fJ2yieldsq0,'LineWidth',1)
plot(I1s,sJ2xi,'-o',I1,fJ2yieldsq,'LineWidth',1)
xlabel('I1')
ylabel('sJ2xi')
%legend('stress path','initial yield surface')
legend('stress path','yield surface')
set(gca,'FontName','Helvetica','FontSize',16)
%
%
% sim=load('ss_enh_isv.txt');
% elem=(sim(:,1));
% loc_flag=(sim(:,2));
% zeta=(sim(:,3));
% gamma_delta=(sim(:,4));
% Q_S=(sim(:,5));
% P_S=(sim(:,6));
% q_St=(sim(:,7));
% cohesion=(sim(:,8));
% friction=(sim(:,9));
% dilation=(sim(:,10));
%
% figure(2)
% plot(zeta,cohesion,'LineWidth',1)
% xlabel('ZETA (m)')
% ylabel('COHESION (MPa) ')
% axis([0 0.005 0 35])
% set(gca,'FontName','Helvetica','FontSize',16)
% %
% figure(3)
% plot(zeta,friction,'LineWidth',1)
% xlabel('ZETA (m)')
% ylabel('FRICTION ANGLE (radian)')
% axis([0 0.005 0 0.55])
% set(gca,'FontName','Helvetica','FontSize',16)
% %
% figure(4)
% plot(zeta,dilation,'LineWidth',1)
% xlabel('ZETA (m)')
% ylabel('DILATION ANGLE (radian)')
% set(gca,'FontName','Helvetica','FontSize',16)
%


