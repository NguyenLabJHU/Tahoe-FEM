
%close all
clear all

format short e

% these stresses are in the intermediate configuration
%  Extract variable ip1.s11 (y/n) ? y
%  Extract variable ip1.s22 (y/n) ? y
%  Extract variable ip1.s33 (y/n) ? y
%  Extract variable ip1.s23 (y/n) ? y
%  Extract variable ip1.s13 (y/n) ? y
%  Extract variable ip1.s12 (y/n) ? y
%  Extract variable ip1.e11 (y/n) ? y
%  Extract variable ip1.e22 (y/n) ? y
%  Extract variable ip1.e33 (y/n) ? y
%  Extract variable ip1.e23 (y/n) ? y
%  Extract variable ip1.e13 (y/n) ? y
%  Extract variable ip1.e12 (y/n) ? y
%  Extract variable ip1.kappa (y/n) ? y
%  Extract variable ip1.plstr (y/n) ? y
%  Extract variable ip1.VM (y/n) ? y
%  Extract variable ip1.press (y/n) ? y
%  Extract variable ip1.ip_loccheck (y/n) ? y
%  Extract variable ip1.el_locflag (y/n) ? y

% Enter GP number "n" 
n=10;
numvars_stress=12;
numvars_isv=6;
num_stress=(n-1)*numvars_stress+2;
num_isv=1+27*numvars_stress+(n-1)*numvars_isv;

sim=load('eldata_DP_SS_361.txt');

time=sim(:,1);

%ip
s11=sim(:,num_stress);
s22=sim(:,num_stress+1);
s33=sim(:,num_stress+2);
s23=sim(:,num_stress+3);
s13=sim(:,num_stress+4);
s12=sim(:,num_stress+5);
e11=sim(:,num_stress+6);
e22=sim(:,num_stress+7);
e33=sim(:,num_stress+8);
e23=sim(:,num_stress+9);
e13=sim(:,num_stress+10);
e12=sim(:,num_stress+11);

kappa=sim(:,num_isv+1);
%final
kappa1=sim(length(time),num_isv+1)
%initial
%kappa0=sim(1,num_isv+1)
kappa0=5.85e-5
%mid
mid_step=round(length(time)/2);
kappa_mid=sim(mid_step,num_isv+1)

plstr=sim(:,num_isv+2);
VM=sim(:,num_isv+3);
press=sim(:,num_isv+4);

beta=0.818;
b=0.818;

Im=[1 0 0; 0 1 0; 0 0 1];
% Principal Stresses
for i=1:length(time)
    sig = [s11(i) s12(i) s13(i) ;
           s12(i) s22(i) s23(i) ;
           s13(i) s23(i) s33(i) ];
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
           e12(i) e22(i) e23(i) ;
           e13(i) e23(i) e33(i) ];
    lambda = eig(eps);
    e1(i) = max(lambda);
    e2(i) = lambda(2);
    e3(i) = min(lambda);
end


pr=linspace(1e-4,-3e-3);
VM0=sqrt(3)*(0.5)*((kappa0-beta*pr)+abs(kappa0-beta*pr));
VM1=sqrt(3)*(0.5)*((kappa0-beta*pr)+abs(kappa0-beta*pr))-kappa1;
VMmid=sqrt(3)*(0.5)*((kappa0-beta*pr)+abs(kappa0-beta*pr))-kappa_mid;


figure(1)
%plot(pr,VM0,'--k',pr,VM1,'-k',press,VM,'-or')
plot(pr,VM0,'--k',pr,VM1,'-k',press,VM,'+b',p,sqrt(3/2)*devSinv,'or')
legend('F0','F','stress path','stress path - calc')
xlabel('mean stress')
ylabel('von Mises stress')
set(gca,'FontName','Helvetica','FontSize',16)

% figure(2)
% plot(pr,devs0inv,'--k',pr,devs1inv,'-k',p,devSinv,'-or')
% legend('F0','Fmid','stress path')
% xlabel('mean stress')
% ylabel('||devS||')
% set(gca,'FontName','Helvetica','FontSize',16)

% figure(3)
% plot(-e33,-s33,'-or','linewidth',1)
% xlabel('axial strain')
% ylabel('axial stress (Pa)')
% set(gca,'FontName','Helvetica','FontSize',16)


