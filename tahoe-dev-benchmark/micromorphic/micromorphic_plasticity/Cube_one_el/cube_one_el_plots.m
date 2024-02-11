
close all
clear all

format short e
% Enter GP number "n" 
n=27;
num=(n-1)*38+2;
% 38+36 (=F,Fe,X,Xe)= 74
%num=(n-1)*74+2;
% 


sim=load('el_Com_Pl_1.txt');
%sim=load('el_cube2x2x2_punch_1.txt');
time=sim(:,1);
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


coh=sim(:,num+18);
hc=sim(:,num+19);
coh_x=sim(:,num+20);
hc_chi=sim(:,num+21);
Delgamma=sim(:,num+22);
Delgammachi=sim(:,num+23);


trsigma=sim(:,num+24);
devsigma=sim(:,num+25);
trrel=sim(:,num+26);
devrel=sim(:,num+27);
trm=sim(:,num+28); %||trm||
devm=sim(:,num+29);

trSPK=sim(:,num+30);
devSPK=sim(:,num+31);
trREL=sim(:,num+32);
devREL=sim(:,num+33);
trM=sim(:,num+34); %||trm||
devM=sim(:,num+35);
invPhi=sim(:,num+36);
invGPhi=sim(:,num+37);


%final
% c=(sim(length(time),20));
 c=(sim(length(time),num+18));
%c=sim(:,num+18); 
% hcf=(sim(length(time),21));
% cx=(sim(length(time),22));
 cx=(sim(length(time),num+20));
%hcfx=(sim(length(time),23));
%  Delgammaf=(sim(length(time),24));
%  Delgammafx=(sim(length(time),25));

%initial
%c0=(sim(1,20));
c0=(sim(1,num+18));
%hc0=(sim(1,21));
%c0x=(sim(1,22));
c0x=(sim(1,num+20));
%hc0x=(sim(1,23));

% Delgamma0=(sim(1,24));
% mid_step=round(length(time)/2);
% cmid=(sim(mid_step,20));
% cmidx=(sim(mid_step,22));
% %
% Delgamma0x=(sim(1,25));
% 
%  for i=1:length(time)
%      cx(i)=1e4;
%      c0x(i)=1e4
%  end

%c0x=1e4;


phi0=0.2;
psi0=0.1;

phi0_x=0.0;
psi0_x=0.0;

beta=-1; %TCs

Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0));
%Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0));
Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0));
Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0));
Bphi=Bphi0;

Aphix=2*sqrt(6)*cx*cos(phi0_x)/(3+beta*sin(phi0_x));
%Aphimidx=2*sqrt(6)*cmidx*cos(phi0_x)/(3+beta*sin(phi0_x));
Aphi0x=2*sqrt(6)*c0x*cos(phi0_x)/(3+beta*sin(phi0_x));
Bphi0x=2*sqrt(6)*sin(phi0_x)/(3+beta*sin(phi0_x));
Bphix=Bphi0x;

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
    %devSinv(i)=norm(devs);
    %devSinv(i)=sqrt(Temp_inv);
    
    prs(i)=trSPK(i)/3;
    prsx(i)=trREL(i)/3;
    
    Temp_inv1=devSPK(i)^2;
    Temp_inv2=devREL(i)^2;
    Stress(i)=sqrt(Temp_inv1+Temp_inv2); 
    %Stress(i)=devSPK(i); 
    %Stress(i)=devREL(i); 
   
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



 pr=linspace(0,-2.5e6);
 prx=linspace(0,-2.5e6);
%  devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
%  devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
 YieldSurface0=(Aphi0-Bphi0*pr+Aphi0x-Bphi0x*prx);
 YieldSurfacefinal=(Aphi-Bphi*pr+Aphix-Bphix*prx);
 %YieldSurface0=(Aphi0-Bphi0*pr+Aphi0x);
 %YieldSurfacefinal=(Aphi-Bphi*pr);

%   YieldSurface0=abs(Aphi0-Bphi0*pr);
%   YieldSurfacefinal=Aphi-Bphi*pr;

% 
% hold on
 figure(4)
plot(Bphi0*pr+Bphi0x*prx,YieldSurface0,'-.g',Bphi*pr+Bphix*prx,YieldSurfacefinal,'-g',Bphi*prs+Bphix*prsx,Stress,'-og');
%plot(Bphi0*pr,YieldSurface0,'--b',Bphi*pr,YieldSurfacefinal,'-.b',Bphi*prs,Stress,'-b');
%plot(pr,YieldSurface0,'--b',pr,YieldSurfacefinal,'-.b',prs,Stress,'-ob');
% legend('F0','F','||devS:devS+devR:devR||')
% xlabel('Bphi*prs+Bphix*prsx')
% ylabel('||devS:devS+devR:devR||')
% set(gca,'FontName','Helvetica','FontSize',16)

% hold on
%  figure(1)
%  plot(time,-s33,'-r','LineWidth',2)
%  xlabel('time','FontSize',18)
%  ylabel('||devM||','FontSize',18)
% h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
%     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
%set(gca,'FontName','Helvetica','FontSize',32)






% 
% figure(1)
% plot(time,coh_x,'-or','LineWidth',2)
% xlabel('time','FontSize',18)
% %ylabel('devrel','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)

% 
% figure(2)
% plot(prsx,devREL,'-or','LineWidth',2)
% xlabel('time','FontSize',18)
% %ylabel('devrel','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)
% % 
% figure(1)
% plot(time,devSPK,'-r','LineWidth',2)
% xlabel('time','FontSize',18)
% ylabel('devSPK','FontSize',18)
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% h=legend('Standard pl.','Macro pl.','Combined pl.','Location','SouthEast')
% set(gca,'FontName','Helvetica','FontSize',32)

% figure(1)
% plot(prsx,devREL,'-r','LineWidth',2)
% xlabel('prsx','FontSize',18)
% ylabel('devREL','FontSize',18)
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)
% 
% 
% figure(4)
% plot(pr,devs0inv,'--k',pr,devsinv,'-k',p,devSinv,'-r');
% legend('F0','F','stress path')
% xlabel('mean stress')
% ylabel('||devS||')
% set(gca,'FontName','Helvetica','FontSize',16)
% 
% %plot(Bphi0*pr+Bphi0x*prx,StressNorm,'--k',pr,devsinv,'-k',p,devSinv,'-r');
% 
% figure(4)
% %plot(Bphi0*pr+Bphi0x*prx,YieldSurface0,'--k',Bphi*pr+Bphix*prx,YieldSurfacefinal,'-k',Bphi*prs+Bphix*prsx,Stress,'-r');
% plot(Bphi*pr+Bphix*prx,YieldSurface0,'--b',Bphi*pr+Bphix*prx,YieldSurfacefinal,'-.b',Bphi*prs+Bphix*prsx,Stress,'-b');
% legend('F0','F','||devS:devS+devR:devR||')
% xlabel('Bphi*pr+Bphix*prx')
% ylabel('||devS:devS+devR:devR||')
% set(gca,'FontName','Helvetica','FontSize',16)



% hold on
% figure(5)
% plot(Bphi*pr+Bphix*prx,YieldSurfacefinal,'-.b',Bphi*prs+Bphix*prsx,Stress,'-b');
% %plot(Bphi0*pr,YieldSurface0,'--b',Bphi*pr,YieldSurfacefinal,'-.b',Bphi*prs,Stress,'-b');
% legend('F','||devS:devS+devR:devR||')
% xlabel('Bphi*prs+Bphix*prsx')
% ylabel('||devS:devS+devR:devR||')
% set(gca,'FontName','Helvetica','FontSize',16)


% figure(4)
% plot(pr,YieldSurface0,'--k',pr,YieldSurfacefinal,'-k',prs,Stress,'-or');
% legend('F0','F','stress path')
% xlabel('pr')
% ylabel('Stress')
% set(gca,'FontName','Helvetica','FontSize',16)

% figure(4)
% plot(Bphi*prs+Bphix*prsx,Stress,'-or');
% legend('F0','F','stress path')
% xlabel('Bphi0*pr+Bphi0x*prx')
% ylabel('Stress')
% set(gca,'FontName','Helvetica','FontSize',16)
% 
% 
% 
% figure(4)
% plot(-Et33,-SPK33,'-or','LineWidth',2)
% xlabel('-FS33','FontSize',18)
% ylabel('-SPK33','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)
% % 
% % 
% figure(5)
% plot(SS33,-REL33,'-or','LineWidth',2)
% xlabel('-SS33','FontSize',18)
% ylabel('-REL33','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)

% figure(2)
% plot(time,Delgamma,'-k','LineWidth',2)
% xlabel('time','FontSize',18)
% 
% %ylabel('cohesion "c"','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)
