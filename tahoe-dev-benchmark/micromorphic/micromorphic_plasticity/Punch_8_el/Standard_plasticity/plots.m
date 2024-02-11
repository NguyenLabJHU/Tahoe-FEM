
%close all
clear all

format short e
% Enter GP number "n" 
n=1;
%num=(n-1)*38+2;
% 38+36 (=F,Fe,X,Xe)= 74
num=(n-1)*74+2;
% 
lambda1=29e6;
mu=7e6;
eta=60e6;
kappa=0.0;
nu=0.0;
sigma_const=0.0;
tau=0.0;

% lambda1=29e6;
% mu=7e6;
% eta=0.0;
% kappa=0.0;
% nu=0.0;
% sigma_const=0.0;
% tau=0.0;
% SPK=zeros(3,3);
Temp_matrix1=zeros(3,3);
Temp_matrix2=zeros(3,3);


%sim=load('el_Col_Pl_1.txt');
sim=load('el_cube2x2x2_punch_7.txt');

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
F11=sim(:,num+38);
F12=sim(:,num+39);
F13=sim(:,num+40);
F21=sim(:,num+41);
F22=sim(:,num+42);
F23=sim(:,num+43);
F31=sim(:,num+44);
F32=sim(:,num+45);
F33=sim(:,num+46);
Fe11=sim(:,num+47);
Fe12=sim(:,num+48);
Fe13=sim(:,num+49);
Fe21=sim(:,num+50);
Fe22=sim(:,num+51);
Fe23=sim(:,num+52);
Fe31=sim(:,num+53);
Fe32=sim(:,num+54);
Fe33=sim(:,num+55);
X11=sim(:,num+56);
X12=sim(:,num+57);
X13=sim(:,num+58);
X21=sim(:,num+59);
X22=sim(:,num+60);
X23=sim(:,num+61);
X31=sim(:,num+62);
X32=sim(:,num+63);
X33=sim(:,num+64);
Xe11=sim(:,num+65);
Xe12=sim(:,num+66);
Xe13=sim(:,num+67);
Xe21=sim(:,num+68);
Xe22=sim(:,num+69);
Xe23=sim(:,num+70);
Xe31=sim(:,num+71);
Xe32=sim(:,num+72);
Xe33=sim(:,num+73);







%final
 c=(sim(length(time),20));
 hcf=(sim(length(time),21));
 cx=(sim(length(time),22));
 hcfx=(sim(length(time),23));
 Delgammaf=(sim(length(time),24));
 Delgammafx=(sim(length(time),25));

%initial
c0=(sim(1,20));
hc0=(sim(1,21));
c0x=(sim(1,22));
hc0x=(sim(1,23));

Delgamma0=(sim(1,24));
mid_step=round(length(time)/2);
cmid=(sim(mid_step,20));
cmidx=(sim(mid_step,22));
%
Delgamma0x=(sim(1,25));










% phi0=30*pi/180;
% psi0=22*pi/180;
% 
% phi0=0.523;
% psi0=0.383;

phi0=0.15;
psi0=0.1;

phi0_x=0.0;
psi0_x=0.0;

beta=-1; %TCs

Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0));
Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0));
Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0));
Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0));
Bphi=Bphi0;

Aphix=2*sqrt(6)*cx*cos(phi0_x)/(3+beta*sin(phi0_x));
Aphimidx=2*sqrt(6)*cmidx*cos(phi0_x)/(3+beta*sin(phi0_x));
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
    devSinv(i)=sqrt(Temp_inv);
    
    prs(i)=trSPK(i)/3;
    prsx(i)=trREL(i)/3;
    
    Temp_inv1=devSPK(i)^2;
    Temp_inv2=devREL(i)^2;
    %Stress(i)=sqrt(Temp_inv1+Temp_inv2); 
    Stress(i)=devSPK(i); 

    
    
    F= [ F11(i) F12(i) F13(i);
         F21(i) F22(i) F23(i);
         F31(i) F32(i) F33(i)];
    Fe= [ Fe11(i) Fe12(i) Fe13(i);
          Fe21(i) Fe22(i) Fe23(i);
          Fe31(i) Fe32(i) Fe33(i)];
   
    Fp=inv(Fe)*F;
    
    if i==1
        Fn_1= Im;
        Fpn_1= Im;
        Xpn_1= Im;   
        Xdn_1=Im;
    else
         Fn_1= [ F11(i-1) F12(i-1) F13(i-1);
              F21(i-1) F22(i-1) F23(i-1);
              F31(i-1) F32(i-1) F33(i-1)];
        
         Fen_1=[ Fe11(i-1) Fe12(i-1) Fe13(i-1);
                 Fe21(i-1) Fe22(i-1) Fe23(i-1);
                 Fe31(i-1) Fe32(i-1) Fe33(i-1)];  
         
         Fpn_1=inv(Fen_1)*Fn_1;
         
         Xdn_1= [ X11(i-1) X12(i-1) X13(i-1);
                  X21(i-1) X22(i-1) X23(i-1);
                  X31(i-1) X32(i-1) X33(i-1)];
         Xen_1= [ Xe11(i-1) Xe12(i-1) Xe13(i-1);
                  Xe21(i-1) Xe22(i-1) Xe23(i-1);
                  Xe31(i-1) Xe32(i-1) Xe33(i-1)];
      
         Xpn_1=inv(Xen_1)*Xdn_1;
                 
    end
   
    Xd= [ X11(i) X12(i) X13(i);
        X21(i) X22(i) X23(i);
        X31(i) X32(i) X33(i)];
    Xe= [ Xe11(i) Xe12(i) Xe13(i);
          Xe21(i) Xe22(i) Xe23(i);
          Xe31(i) Xe32(i) Xe33(i)];
      
   Xp=inv(Xe)*Xd;
  
   deltat=time(2)-time(1);
   
   Ce=transpose(Fe)*Fe;
   Ct=transpose(F)*F;
   Cxe=transpose(Xe)*Xe;
   if i==1
        Lp=(1/deltat)*(Fp)*inv(Fp);
       Lxp=(1/deltat)*(Xp)*inv(Xp);
   else
        Lp=(1/deltat)*(Fp-Fpn_1)*inv(Fp);
       Lxp=(1/deltat)*(Xp-Xpn_1)*inv(Xp); 
   end
   
   
   Et=0.5*(Ct-Im);

   CeLp=Ce*Lp;
   PSIe=transpose(Fe)*Xe;
   PSIt=transpose(F)*Xd;
   
   Ee=0.5*(Ce-Im);
   trEe=Ee(1,1)+Ee(2,2)+Ee(3,3);
   Epsilone=PSIe-Im;
   Epsilont=PSIt-Im;
   
   trEp=Epsilone(1,1)+Epsilone(2,2)+Epsilone(3,3);
   SST=PSIe*Lxp*inv(Cxe)*transpose(PSIe);
   
   SPK=Im*trEe*(lambda1+tau)+Ee*2*(sigma_const+mu)+eta*trEp*Im+kappa*Epsilone+nu*transpose(Epsilone);
   REL=tau*trEe*Im+2*sigma_const*Ee+(eta-tau)*trEp*Im+(nu-sigma_const)*Epsilone+(kappa-sigma_const)*transpose(Epsilone);
   
    SPK33(i)=SPK(3,3);
    REL33(i)=REL(3,3);
    Et33(i)=Et(3,3);
    Epsilont33(i)=Epsilont(3,3);

   
   % S(KL) [F(kK)Fdot(kL)]= S(KL) [F(kK)(1/deltat)(F-Fn)(kL)] 
   %(SIGMA-S)(KL) [F(kK)Xdot(kA)X^(-1)(Al)F(lL)]=(SIGMA-S)(KL) [F(kK)(1/deltat)(X-Xn)(kA)X^(-1)(Al)F(lL)]
   
   FTFdot=(1/deltat)*F'*(F-Fn_1);
   FS33(i)=FTFdot(3,3);
   
   FTXdotX1F=(1/deltat)*F'*(Xd-Xdn_1)*inv(Xd)*F;
   SS33(i)=FTXdotX1F(3,3);
   
    
    SPK33(i)=SPK(3,3);
    REL33(i)=REL(3,3);
%    
%    CeLp33(i)=CeLp(3,3);
%    SST33(i)=SST(3,3);
   
    
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



 pr=linspace(0,-10e5);
 prx=linspace(0,-10e5);
%  devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
%  devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
%  YieldSurface0=Aphi0-Bphi0*pr+Aphi0x-Bphi0x*prx;
%  YieldSurfacefinal=Aphi-Bphi*pr+Aphix-Bphix*prx;
  YieldSurface0=abs(Aphi0-Bphi0*pr);
  YieldSurfacefinal=Aphi-Bphi*pr;

% 
% figure(1)
% plot(time,coh,'-or','LineWidth',2)
% xlabel('time','FontSize',18)
% %ylabel('devrel','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)
% 
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
% 
% figure(3)
% plot(prs,devSPK,'-or','LineWidth',2)
% xlabel('prs','FontSize',18)
% %ylabel('devrel','FontSize',18)
% % h=legend('Micromorphic FSE with \lambda, \mu, \kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu','Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma'...
% %     ,'Micromorphic FSE with \lambda, \mu, \eta,\kappa,\nu,\tau,\sigma and \tau_{7 }','Location','NorthWest');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',32)


% figure(4)
% plot(pr,devs0inv,'--k',pr,devsinv,'-k',p,devSinv,'-r');
% legend('F0','F','stress path')
% xlabel('mean stress')
% ylabel('||devS||')
% set(gca,'FontName','Helvetica','FontSize',16)

%plot(Bphi0*pr+Bphi0x*prx,StressNorm,'--k',pr,devsinv,'-k',p,devSinv,'-r');
% 
%hold on
figure(4)
plot(pr,YieldSurface0,'-xr',pr,YieldSurfacefinal,'-+k',prs,Stress,'-or');
legend('F0','F','||devS||')
xlabel('pr')
ylabel('||devS||')
set(gca,'FontName','Helvetica','FontSize',16)
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
