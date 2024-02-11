
%close all
clear all

format short e

sim=load('eldataout_1.txt');
time=sim(:,1);
s11=(sim(:,2))/1000;
s22=(sim(:,3))/1000;
s33=(sim(:,4))/1000;
s23=(sim(:,5))/1000;
s13=(sim(:,6))/1000;
s12=(sim(:,7))/1000;
e11=sim(:,8);
e22=sim(:,9);
e33=sim(:,10);
e23=sim(:,11);
e13=sim(:,12);
e12=sim(:,13);
%final
kappa=(sim(length(time),14))/1000;
c=(sim(length(time),15))/1000;
%initial
kappa0=(sim(1,14))/1000;
c0=(sim(1,15))/1000;
%mid
mid_step=round(length(time)*(2/3));
kappamid=(sim(mid_step,14))/1000;
cmid=(sim(mid_step,15))/1000;
%
p_prime=(sim(:,16))/1000;
I1_prime=3*p_prime;
sdev_sdev=(sim(:,17))/(1e6);
sqrt_sdev_sdev = sqrt(sdev_sdev);

phi0=pi*20/180; %radians
psi0=pi*10/180; %radians

R=10; %dimensionless

beta=-1; %TCs

Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0))
Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0))
Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0))
Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0))

chi0=10;


% Principal Stresses
for i=1:length(time)
    sig = [s11(i) s12(i) s13(i) ;
           s12(i) s22(i) s23(i) ;
           s13(i) s23(i) s33(i) ];
    lambda = eig(sig);
    s1(i,1) = max(lambda);
    s2(i,1) = lambda(2);
    s3(i,1) = min(lambda);
end

I1=s1+s2+s3;
meanstress=I1./3;

% Principal Strains
for i=1:length(time)
    eps = [e11(i) e12(i) e13(i) ;
           e12(i) e22(i) e23(i) ;
           e13(i) e23(i) e33(i) ];
    lambda = eig(eps);
    e1(i,1) = max(lambda);
    e2(i,1) = lambda(2);
    e3(i,1) = min(lambda);
end



%shear yield surface (initial, and final)
X0=kappa0-R*(Aphi0-Bphi0*kappa0)
X=kappa-R*(Aphi-Bphi0*kappa)
I1f=X:.1:chi0;
I1f0=X0:.1:chi0;
FyI0=Aphi0-Bphi0*I1f0/3;
FyI0=0.5*(abs(FyI0)+FyI0);
FyI=Aphi-Bphi0*I1f/3;
FyI=0.5*(abs(FyI)+FyI);
%cap surface
func10=(abs(kappa0-I1f0)+(kappa0-I1f0))/2;
Fcf0=1-func10.*(kappa0-I1f0)/((X0-kappa0)^2);
func1=(abs(kappa-I1f)+(kappa-I1f))/2;
Fcf=1-func1.*(kappa-I1f)/((X-kappa)^2);
%Fcf=0.5*(abs(Fcf)+Fcf);
J20=0.5*Fcf0.*FyI0.^2;
ssf0=2*J20;
snormf0=sqrt(ssf0);
J2sqf0=sqrt(J20);
J2=0.5*Fcf.*FyI.^2;
ssf=2*J2;
snormf=sqrt(ssf);
J2sqf=sqrt(J2);

%plastic potential
Apsi=2*sqrt(6)*c*cos(psi0)/(3+beta*sin(psi0))
Bpsi=2*sqrt(6)*sin(psi0)/(3+beta*sin(psi0))

%shear plastic potential surface
Xg=kappa-R*(Apsi-Bpsi*kappa)
I1g=Xg:.1:chi0;
FyI=Apsi-Bpsi*I1g/3;
FyI=0.5*(abs(FyI)+FyI);
%cap surface
func1=(abs(kappa-I1g)+(kappa-I1g))/2;
Fcg=1-func1.*(kappa-I1g)/((Xg-kappa)^2);
Fcg=0.5*(abs(Fcg)+Fcg);
J2=0.5*Fcg.*FyI.^2;
ssg=2.0*J2;
snormg=sqrt(ssg);
J2sqg=sqrt(J2);

%mid step
Xmid=kappamid-R*(Aphimid-Bphi0*kappamid)
I1fmid=Xmid:.1:chi0;
FyI=Aphimid-Bphi0*I1fmid/3;
FyI=0.5*(abs(FyI)+FyI);
%cap surface
func1=(abs(kappamid-I1fmid)+(kappamid-I1fmid))/2;
Fcf=1-func1.*(kappamid-I1fmid)/((Xmid-kappamid)^2);
J2=0.5*Fcf.*FyI.^2;
ssf=2*J2;
snormfmid=sqrt(ssf);
J2sqf=sqrt(J2);


% stress strain
figure(1)
plot(abs(e33),abs(s33-s11),'-k','LineWidth',2)
%plot(abs(e33),abs(s33),'-k','LineWidth',2)
%plot(abs(e11),abs(s11),'--k',abs(e33),abs(s33),'-k','LineWidth',2)
xlabel('e33')
ylabel('s33-s11')
%ylabel('s33')
%ylabel('s')
%legend('s11-e11','s33-e33')
set(gca,'FontName','Helvetica','FontSize',16)
set(gca,'XTickLabel',{'0','1','2','3','4','5'})

% principal stress strain
% figure(3)
% %plot(abs(e3),abs(s3),'LineWidth',2)
% plot(abs(e3),abs(s3),'-k','LineWidth',2)
% xlabel('e3')
% ylabel('s3')
% set(gca,'FontName','Helvetica','FontSize',16)


% figure(2)
% plot(I1f,ssf,'LineWidth',2)
% xlabel('I1')
% ylabel('2J2')
% %legend('yield surface')
% set(gca,'FontName','Helvetica','FontSize',16)



% comparison yield and plastic potential surface
figure(2)
%plot(I1g,J2sqg,I1f,J2sqf)
%plot(I1f,snormf,'-k',I1fmid,snormfmid,':k',I1f0,snormf0,'--k',I1g,snormg,'-.k',I1_prime,sqrt_sdev_sdev,'ok','LineWidth',2)
plot(I1f,snormf,'-k',I1fmid,snormfmid,':k',I1f0,snormf0,'--k',I1_prime,sqrt_sdev_sdev,'-.ok','LineWidth',2)
%plot(I1f,snormf,'-k',I1f0,snormf0,'--k',I1g,snormg,'-.k',I1_prime,sqrt_sdev_sdev,'ok','LineWidth',2)
%plot(I1f,snormf,'-k',I1f0,snormf0,'--k',I1_prime,sqrt_sdev_sdev,'ok','LineWidth',2)
xlabel('I1')
ylabel('s2J2')
%legend('Ff','Fmid','F0','Gf','sim')
legend('F','Fmid','F0','stress path')
axis([-200 0 0 20])
%legend('Ff','Fi','Gf','sim')
%legend('Ff','Fi','sim')
set(gca,'FontName','Helvetica','FontSize',16)

%
% figure(2)
% %plot(I1g,J2sqg,I1f,J2sqf)
% plot(I1f,ssf,'-k',I1f,ssf0,'--k',I1g,ssg,'-.k',I1_prime,sdev_sdev,'ok','LineWidth',2)
% xlabel('I1')
% ylabel('2J2')
% legend('F','Fi','G','sim')
% set(gca,'FontName','Helvetica','FontSize',16)

% cap function
% figure(3)
% plot(I1f,Fcf)
% legend('Fcapphi')
% grid on

% % volumetric plastic flow
% Cpsi = -(2*func1.*(Apsi - Bpsi*I1g/3).^2)/((X - kappa0)^2) + (2*Bpsi/3)*Fc.*(Apsi - Bpsi*I1g/3);
% figure(3)
% plot(I1g,Cpsi)
% legend('Cpsi')
% grid on
% 
% %hardening function for kappa
% dhdkappa = 2*func1.*(R*Bpsi*(kappa0-I1g)/(X-kappa0) - 1).*(FyI.^2)/((X-kappa0)^2);
% figure(4)
% plot(I1g,dhdkappa)
% legend('dhdkappa')
% grid on
% 
% %hardening function for cohesion
% dhdc = 2*(Apsi_noc)*Fc.*FyI;
% figure(5)
% plot(I1g,dhdc)
% legend('dhdc')
% grid on
% 
% alternate evolution function for kappa
% %epsvol=-1:.01:0.1;
% epsvol=0:.001:1.0;
% signfunc=1.0;
% alphakap=100;
% %kappa = signfunc*kappa0*exp(-alphakap*epsvol);
% kappafunc = signfunc*(-kappa0)*exp(alphakap*epsvol);
% figure(5)
% plot(epsvol,kappafunc,'-k','LineWidth',2)
% axis([0 .1 0 5e4])
% xlabel('epsvol')
% ylabel('kappa')
% %legend('kappa')
% set(gca,'FontName','Helvetica','FontSize',16)
% set(gca,'YTickLabel',{'0','10','20','30','40','50'})
% %grid on

