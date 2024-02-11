
close all
clear all

format short e

% Enter GP number "n" 
n=1;
numvars=22;
num=(n-1)*numvars+2;

sim=load('eldata_cap_1.txt');

time=sim(:,1);

%  Extract variable ip1.s11 (y/n) ? y
%  Extract variable ip1.s22 (y/n) ? y
%  Extract variable ip1.s33 (y/n) ? y
%  Extract variable ip1.s23 (y/n) ? y
%  Extract variable ip1.s13 (y/n) ? y
%  Extract variable ip1.s12 (y/n) ? y
%  Extract variable ip1.p_f (y/n) ? y
%  Extract variable ip1.e11 (y/n) ? y
%  Extract variable ip1.e22 (y/n) ? y
%  Extract variable ip1.e33 (y/n) ? y
%  Extract variable ip1.e23 (y/n) ? y
%  Extract variable ip1.e13 (y/n) ? y
%  Extract variable ip1.e12 (y/n) ? y
%  Extract variable ip1.phi_s (y/n) ? y
%  Extract variable ip1.phi_f (y/n) ? y
%  Extract variable ip1.J (y/n) ? y
%  Extract variable ip1.k (y/n) ? y
%  Extract variable ip1.kappa (y/n) ? y
%  Extract variable ip1.c (y/n) ? y
%  Extract variable ip1.p_prime (y/n) ? y
%  Extract variable ip1.sdev_sdev (y/n) ? y
%  Extract variable ip1.eps_vol_p (y/n) ? y

%ip
s11=sim(:,num);
s22=sim(:,num+1);
s33=sim(:,num+2);
s23=sim(:,num+3);
s13=sim(:,num+4);
s12=sim(:,num+5);
p_f=sim(:,num+6);
e11=sim(:,num+7);
e22=sim(:,num+8);
e33=sim(:,num+9);
e23=sim(:,num+10);
e13=sim(:,num+11);
e12=sim(:,num+12);

phi_s=sim(:,num+13);
phi_f=sim(:,num+14);
Jac=sim(:,num+15);
k_perm=sim(:,num+16);

%final
kappa=sim(length(time),num+17)
c=sim(length(time),num+18)
%initial
kappa0=sim(1,num+17)
c0=sim(1,num+18)
%mid
mid_step=round(length(time)*(2/3));
kappamid=sim(mid_step,num+17);
cmid=sim(mid_step,num+18);

p_prime=sim(:,num+19);
I1_prime=3*p_prime;
sdev_sdev=sim(:,num+20);
sqrt_sdev_sdev = sqrt(sdev_sdev);

eps_vol_p=sim(:,num+21);

phi0=0.61; %radians
psi0=0.61; %radians

R=10; %dimensionless

beta=-1; %TCs

Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0))
Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0))
Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0))
Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0))

chi0=1e5;
chi_inc=1e3;

Im=[1 0 0; 0 1 0; 0 0 1];
% Principal Stresses
for i=1:length(time)
    sig = [s11(i) s12(i) s13(i) ;
           s12(i) s22(i) s23(i) ;
           s13(i) s23(i) s33(i) ];
    p(i)=(sig(1,1)+sig(2,2)+sig(3,3))/3;   
    devs=sig-p(i)*Im;
    devSinv_norm(i)=normest(devs,2);
    Temp_inv1=devs(1,1)*devs(1,1)+devs(2,2)*devs(2,2)+devs(3,3)*devs(3,3)+...
        2*(devs(1,3)*devs(1,3)+devs(1,2)*devs(1,2)+devs(2,3)*devs(2,3));
    devSinv1(i)=sqrt(Temp_inv1);
    Temp_inv2=0;
    for j=1:3
        for k=1:3
            Temp_inv2=Temp_inv2+devs(j,k)*devs(j,k);
        end
    end 
    devSinv2(i)=sqrt(Temp_inv2);
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
I1f=X:chi_inc:chi0;
I1f0=X0:chi_inc:chi0;
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
I1g=Xg:chi_inc:chi0;
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
I1fmid=Xmid:chi_inc:chi0;
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
% figure(1)
% %plot(abs(e33),abs(s33-s11),'-k','LineWidth',2)
% %plot(-(e22),-(s22),'-k','LineWidth',2)
% plot(-(e33),-(s33),'-k','LineWidth',2)
% %plot(e22,s22,'-k','LineWidth',2)
% %plot(abs(e11),abs(s11),'--k',abs(e33),abs(s33),'-k','LineWidth',2)
% xlabel('e33')
% %ylabel('s33-s11')
% ylabel('s33')
% %ylabel('s')
% %legend('s11-e11','s33-e33')
% %set(gca,'XTickLabel',{'0','0.5','1','1.5','2','2.5','3'})
% set(gca,'FontName','Helvetica','FontSize',16)

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
%plot(I1f,snormf,'-k',I1fmid,snormfmid,':k',I1f0,snormf0,'--k',I1_prime,sqrt_sdev_sdev,'.k','LineWidth',1)
%plot(I1f,snormf,'-k',I1fmid,snormfmid,':k',I1f0,snormf0,'--k',I1_prime,sqrt_sdev_sdev,'.k',...
%    p,devSinv_norm,'ok',p,devSinv1,'sk',p,devSinv2,'+k','LineWidth',1)
plot(I1f,snormf,'-k',I1fmid,snormfmid,':k',I1f0,snormf0,'--k',I1_prime,sqrt_sdev_sdev,'.k',...
    p,devSinv_norm,'ok','LineWidth',1)
%plot(I1f,snormf,'-k',I1f0,snormf0,'--k',I1g,snormg,'-.k',I1_prime,sqrt_sdev_sdev,'ok','LineWidth',2)
%plot(I1f,snormf,'-k',I1f0,snormf0,'--k',I1_prime,sqrt_sdev_sdev,'ok','LineWidth',2)
xlabel('I1')
ylabel('s2J2')
%legend('Ff','Fmid','F0','Gf','sim')
%legend('F','Fmid','F0','sim','sim-norm','sim-sqrt1','sim-sqrt2')
legend('F','Fmid','F0','sim - inter config','sim - curr config')
%axis([-100 0 0 12])  
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

