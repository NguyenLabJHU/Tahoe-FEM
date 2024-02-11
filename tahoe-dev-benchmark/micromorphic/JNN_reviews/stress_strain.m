
%close all
clear all

format short e

% these are stresses in the intermediate configuration
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

sim=load('eldata_1.txt');

time=sim(:,1);

% Enter GP number "n" 
n=27;
numvars=35;
num=(n-1)*numvars+2;

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

kc=sim(:,num+18);
khc=sim(:,num+19);
kDelgamma=sim(:,num+20);
trSigma=sim(:,num+21);
devSigma=sim(:,num+22);
trRel=sim(:,num+23);
devRel=sim(:,num+24);
trm=sim(:,num+25);
devm=sim(:,num+26);
trS=sim(:,num+27);
devS=sim(:,num+28);
trRelS=sim(:,num+29);
devRelS=sim(:,num+30);
trM=sim(:,num+31);
devM=sim(:,num+32);
PHI=sim(:,num+33);
GPHI=sim(:,num+34);

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

stressfactor=1e6;

figure(1)
plot(time,s33*stressfactor,'-or',time,s22*stressfactor,'-ob')
legend('s33','s22')
xlabel('time')
ylabel('stress (MPa)')
grid on
set(gca,'FontName','Helvetica','FontSize',16)

figure(2)
%plot(time,devS*stressfactor,'-or',time,devSigma*stressfactor,'-ob',time,devSinv*stressfactor,'-sg')
plot(time,devS*stressfactor,'-or',time,devSigma*stressfactor,'-ob')
legend('||devS||','||devs||')
xlabel('time')
ylabel('stress (MPa)')
grid on
set(gca,'FontName','Helvetica','FontSize',16)

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

% sim_drain=load('eldata_drain_1.txt');
% 
% time_drain=sim_drain(:,1);
% 
% % Enter GP number "n" 
% n=25;
% numvars=22;
% num=(n-1)*numvars+2;
% 
% %ip
% s11_drain=sim_drain(:,num);
% s22_drain=sim_drain(:,num+1);
% s33_drain=sim_drain(:,num+2);
% s23_drain=sim_drain(:,num+3);
% s13_drain=sim_drain(:,num+4);
% s12_drain=sim_drain(:,num+5);
% e11_drain=sim_drain(:,num+7);
% e22_drain=sim_drain(:,num+8);
% e33_drain=sim_drain(:,num+9);
% e23_drain=sim_drain(:,num+10);
% e13_drain=sim_drain(:,num+11);
% e12_drain=sim_drain(:,num+12);
% 
% figure(3)
% plot(time,s33*stressfactor,'-or',time_drain,s33_drain*stressfactor,'-ob')
% legend('s33','s33drain')
% xlabel('time')
% ylabel('stress (MPa)')
% grid on
% set(gca,'FontName','Helvetica','FontSize',16)

