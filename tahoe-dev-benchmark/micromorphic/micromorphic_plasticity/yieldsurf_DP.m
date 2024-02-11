
%close all
clear all

format short e

% these stresses are in the intermediate configuration
%  Extract variable ip1.s11 (y/n) ? y
%  Extract variable ip1.s22 (y/n) ? y
%  Extract variable ip1.s33 (y/n) ? y
%  Extract variable ip1.s12 (y/n) ? y
%  Extract variable ip1.s13 (y/n) ? y
%  Extract variable ip1.s21 (y/n) ? y
%  Extract variable ip1.s23 (y/n) ? y
%  Extract variable ip1.s31 (y/n) ? y
%  Extract variable ip1.s32 (y/n) ? y
%  Extract variable ip1.e11 (y/n) ? y
%  Extract variable ip1.e22 (y/n) ? y
%  Extract variable ip1.e33 (y/n) ? y
%  Extract variable ip1.e12 (y/n) ? y
%  Extract variable ip1.e13 (y/n) ? y
%  Extract variable ip1.e21 (y/n) ? y
%  Extract variable ip1.e23 (y/n) ? y
%  Extract variable ip1.e31 (y/n) ? y
%  Extract variable ip1.e32 (y/n) ? y
%  Extract variable ip1.kc (y/n) ? y
%  Extract variable ip1.khc (y/n) ? y
%  Extract variable ip1.kDelgamma (y/n) ? y
%  Extract variable ip1.trSigma (y/n) ? y
%  Extract variable ip1.||dev(Sigma)|| (y/n) ? y
%  Extract variable ip1.trRel (y/n) ? y
%  Extract variable ip1.||dev(Rel)|| (y/n) ? y
%  Extract variable ip1.trm (y/n) ? y
%  Extract variable ip1.||dev(m)|| (y/n) ? y
%  Extract variable ip1.trS (y/n) ? y
%  Extract variable ip1.||dev(S)|| (y/n) ? y
%  Extract variable ip1.trRelS (y/n) ? y
%  Extract variable ip1.||dev(RelS)|| (y/n) ? y
%  Extract variable ip1.trM (y/n) ? y
%  Extract variable ip1.||dev(M)|| (y/n) ? y
%  Extract variable ip1.||PHI|| (y/n) ? y
%  Extract variable ip1.||GRAD(PHI)|| (y/n) ? y

% Enter GP number "n" 
n=25;
numvars=35
num=(n-1)*numvars+2;

sim=load('eldata_1.txt');

time=sim(:,1);

%ip
s11=sim(:,num);
s22=sim(:,num+1);
s33=sim(:,num+2);
s12=sim(:,num+3);
s13=sim(:,num+4);
s21=sim(:,num+5);
s23=sim(:,num+6);
s31=sim(:,num+7);
s32=sim(:,num+8);
e11=sim(:,num+9);
e22=sim(:,num+10);
e33=sim(:,num+11);
e12=sim(:,num+12);
e13=sim(:,num+13);
e21=sim(:,num+14);
e23=sim(:,num+15);
e31=sim(:,num+16);
e32=sim(:,num+17);

%final
c=sim(length(time),num+18)
hc=sim(length(time),num+19);
Delgamma=sim(length(time),num+20)
%initial
c0=sim(1,num+18)
hc0=sim(1,num+19);
Delgamma0=sim(1,num+20)
%
mid_step=round(length(time)/2);
cmid=sim(mid_step,num+18);

c_all=sim(:,num+18);
Delgamma_all=sim(:,num+20);

phi0=0.61;
psi0=0.61;

beta=-1; %TCs

Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0))
Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0))
Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0))
Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0))

Im=[1 0 0; 0 1 0; 0 0 1];
% Principal Stresses
for i=1:length(time)
    sig = [s11(i) s12(i) s13(i) ;
           s21(i) s22(i) s23(i) ;
           s31(i) s32(i) s33(i) ];
    lambda = eig(sig);
    s1(i) = max(lambda);
    s2(i) = lambda(2);
    s3(i) = min(lambda);
    p(i)=(sig(1,1)+sig(2,2)+sig(3,3))/3;
    devs=sig-p(i)*Im;
    Temp_inv=0;
    for j=1:3
        for k=1:3
            Temp_inv=Temp_inv+devs(j,k)*devs(j,k);
        end
    end 
    %devSinv(i)=normest(devs,2);
    devSinv(i)=sqrt(Temp_inv);
end

I1=s1+s2+s3;
meanstress=I1./3;

% Principal Strains
for i=1:length(time)
    eps = [e11(i) e12(i) e13(i) ;
           e21(i) e22(i) e23(i) ;
           e31(i) e32(i) e33(i) ];
    lambda = eig(eps);
    e1(i) = max(lambda);
    e2(i) = lambda(2);
    e3(i) = min(lambda);
end


pr=linspace(0,-2e6);
devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
devs1inv=0.5*((Aphimid-Bphi0*pr)+abs(Aphimid-Bphi0*pr));
devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));


figure(1)
plot(pr,devs0inv,'--k',pr,devsinv,'-k',p,devSinv,'-or')
legend('F0','F','stress path')
xlabel('mean stress')
ylabel('||devS||')
set(gca,'FontName','Helvetica','FontSize',16)

% figure(2)
% plot(pr,devs0inv,'--k',pr,devs1inv,'-k',p,devSinv,'-or')
% legend('F0','Fmid','stress path')
% xlabel('mean stress')
% ylabel('||devS||')
% set(gca,'FontName','Helvetica','FontSize',16)

figure(3)
plot(-e33,-s33,'-or','linewidth',1)
xlabel('axial strain')
ylabel('axial stress (Pa)')
set(gca,'FontName','Helvetica','FontSize',16)


