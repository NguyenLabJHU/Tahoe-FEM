clear all
%
sim=load('eldataloc_1.txt');
time=(sim(:,1));
sig11=(sim(:,2));
sig22=(sim(:,3));
sig33=(sim(:,4));
sig23=(sim(:,5));
sig13=(sim(:,6));
sig12=(sim(:,7));
eps11=(sim(:,8));
eps22=(sim(:,9));
eps33=(sim(:,10));
eps23=(sim(:,11));
eps13=(sim(:,12));
eps12=(sim(:,13));
loc_ip1=(sim(:,14));
loc_ip2=(sim(:,15));
loc_ip3=(sim(:,16));
loc_ip4=(sim(:,17));
loc_ip5=(sim(:,18));
loc_ip6=(sim(:,19));
loc_ip7=(sim(:,20));
%ip8
alph11=(sim(:,21));
alph22=(sim(:,22));
alph33=(sim(:,23));
alph23=(sim(:,24));
alph13=(sim(:,25));
alph12=(sim(:,26));
kappa=(sim(:,27));
press=sim(:,28);
J2=(sim(:,29));
J3=(sim(:,30));
loc_ip8=(sim(:,31));
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
D1=1.467e-3;
D2=0.0;
W=.08;

X0=-460; % MPa
kappa0=-8.05; %MPa
I1=-500:1:50;
%shear failure surface
Ff=A-C*exp(B*I1)-theta*I1;
%Ff=A-C*exp(B*I1);
Ff=0.5*(abs(Ff)+Ff);
%cap surface
X=kappa0-R*(A-C*exp(B*kappa0)-theta*kappa0);
Fc=1+(I1-kappa0).*(abs(I1-kappa0)-(I1-kappa0))/(2*(X-kappa0)^2);
Fc=0.5*(abs(Fc)+Fc);
Fy=0.5*(abs(Ff-N)+Ff-N);
fJ2xi=Fc.*Fy.^2;
fJ2xisq=sqrt(Fc).*Fy;
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
plot(abs(eps22),abs(sig22))
xlabel('STRAIN')
ylabel('STRESS (MPa)')
%legend('100 steps','40 steps')
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
figure(8)
plot(I1s,sJ2,'o')
xlabel('I1')
ylabel('sJ2')
%
figure(9)
plot(I1s,sJ2xi,'o',I1,fJ2xisq)
xlabel('I1')
ylabel('sJ2xi')
legend('sJ2xi','yield')
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
ylabel('COHESION (MPa) ')
axis([0 0.005 0 35])
%
figure(3)
plot(zeta,friction)
xlabel('ZETA (m)')
ylabel('FRICTION ANGLE (radian)')
axis([0 0.005 0 0.55])
%
figure(4)
plot(zeta,dilation)
xlabel('ZETA (m)')
ylabel('DILATION ANGLE (radian)')
%


