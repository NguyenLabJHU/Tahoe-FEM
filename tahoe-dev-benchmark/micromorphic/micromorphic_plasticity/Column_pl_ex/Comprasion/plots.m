
close all
clear all

format short e
% Enter GP number "n" 
n=1;
%num=(n-1)*38+2;
% 38+36 (=F,Fe,X,Xe)= 74
num=(n-1)*74+2;
% 

sim=load('el_Comb_Pl_1.txt');
timecomb=sim(:,1);

s33comb=sim(:,num+2);
e33comb=sim(:,num+11);


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
devSPKcomb=sim(:,num+31);
trREL=sim(:,num+32);
devREL=sim(:,num+33);
trM=sim(:,num+34); %||trm||
devM=sim(:,num+35);
invPhi=sim(:,num+36);
invGPhi=sim(:,num+37);

%final
 c=(sim(length(timecomb),20));
 hcf=(sim(length(timecomb),21));
 cx=(sim(length(timecomb),22));
 hcfx=(sim(length(timecomb),23));
 Delgammaf=(sim(length(timecomb),24));
 Delgammafx=(sim(length(timecomb),25));

%initial
c0=(sim(1,20));
hc0=(sim(1,21));
c0x=(sim(1,22));
hc0x=(sim(1,23));

Delgamma0=(sim(1,24));
mid_step=round(length(timecomb)/2);
cmid=(sim(mid_step,20));
cmidx=(sim(mid_step,22));
%
Delgamma0x=(sim(1,25));


phi0=0.15;
psi0=0.1;

phi0_x=0.15;
psi0_x=0.05;

beta=-1; %TCs

Aphicomb=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0));
Aphimidcomb=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0));
Aphi0comb=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0));
Bphi0comb=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0));
Bphicomb=Bphi0comb;

Aphixcomb=2*sqrt(6)*cx*cos(phi0_x)/(3+beta*sin(phi0_x));
Aphimidxcomb=2*sqrt(6)*cmidx*cos(phi0_x)/(3+beta*sin(phi0_x));
Aphi0xcomb=2*sqrt(6)*c0x*cos(phi0_x)/(3+beta*sin(phi0_x));
Bphi0xcomb=2*sqrt(6)*sin(phi0_x)/(3+beta*sin(phi0_x));
Bphixcomb=Bphi0xcomb;

% Principal Stresses
for i=1:length(timecomb)
   
    prcomb(i)=trSPK(i)/3;
    prsxcomb(i)=trREL(i)/3;
    
    Temp_inv1=devSPKcomb(i)^2;
    Temp_inv2=devREL(i)^2;
    stresscomb(i)=sqrt(Temp_inv1+Temp_inv2); 
    %Stress(i)=devSPK(i); 

end




 prscomb=linspace(0,-50e4);
 prxcomb=linspace(0,-50e4);
%  devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
%  devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
 YieldSurface0comb=Aphi0comb-Bphi0comb*prscomb+Aphi0xcomb-Bphi0xcomb*prxcomb;
 YieldSurfacecomb=Aphicomb-Bphicomb*prscomb+Aphixcomb-Bphixcomb*prxcomb;
%   YieldSurface0=abs(Aphi0-Bphi0*pr);
%   YieldSurfacefinal=Aphi-Bphi*pr;



sim=load('el_Standard_Pl_1.txt');
times=sim(:,1);

s33standard=sim(:,num+2);
e33standard=sim(:,num+11);


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
devSPKs=sim(:,num+31);
trREL=sim(:,num+32);
devREL=sim(:,num+33);
trM=sim(:,num+34); %||trm||
devM=sim(:,num+35);
invPhi=sim(:,num+36);
invGPhi=sim(:,num+37);

%final
 c=(sim(length(times),20));
 hcf=(sim(length(times),21));
 cx=(sim(length(times),22));
 hcfx=(sim(length(times),23));
 Delgammaf=(sim(length(times),24));
 Delgammafx=(sim(length(times),25));

%initial
c0=(sim(1,20));
hc0=(sim(1,21));
c0x=(sim(1,22));
hc0x=(sim(1,23));

Delgamma0=(sim(1,24));
mid_step=round(length(times)/2);
cmid=(sim(mid_step,20));
cmidx=(sim(mid_step,22));
%
Delgamma0x=(sim(1,25));


phi0=0.15;
psi0=0.1;

phi0_x=0.15;
psi0_x=0.05;

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

% Principal Stresses
for i=1:length(times)
   
    prstandard(i)=trSPK(i)/3;
    %prsx(i)=trREL(i)/3;
    
    Temp_inv1=devSPKs(i)^2;
    %Temp_inv2=devREL(i)^2;
    %Stress(i)=sqrt(Temp_inv1+Temp_inv2); 
    stressstandard(i)=devSPKs(i); 

end




 prsstandard=linspace(0,-100e4);
 prx=linspace(0,-10e4);
%  devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
%  devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
%  YieldSurface0standard=Aphi0-Bphi0*pr+Aphi0x-Bphi0x*prx;
%  w=Aphi-Bphi*pr+Aphix-Bphix*prx;
  YieldSurface0standard=abs(Aphi0-Bphi0*prsstandard);
  YieldSurfacestandard=Aphi-Bphi*prsstandard;

  
sim=load('el_Macro_Pl_1.txt');
timeM=sim(:,1);

s33Macro=sim(:,num+2);
e33Macro=sim(:,num+11);


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
devSPKM=sim(:,num+31);
trREL=sim(:,num+32);
devREL=sim(:,num+33);
trM=sim(:,num+34); %||trm||
devM=sim(:,num+35);
invPhi=sim(:,num+36);
invGPhi=sim(:,num+37);

%final
 c=(sim(length(timeM),20));
 hcf=(sim(length(timeM),21));
 cx=(sim(length(timeM),22));
 hcfx=(sim(length(timeM),23));
 Delgammaf=(sim(length(timeM),24));
 Delgammafx=(sim(length(timeM),25));

%initial
c0=(sim(1,20));
hc0=(sim(1,21));
c0x=(sim(1,22));
hc0x=(sim(1,23));

Delgamma0=(sim(1,24));
mid_step=round(length(timeM)/2);
cmid=(sim(mid_step,20));
cmidx=(sim(mid_step,22));
%
Delgamma0x=(sim(1,25));


phi0=0.15;
psi0=0.1;

phi0_x=0.15;
psi0_x=0.05;

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

% Principal Stresses
for i=1:length(timeM)
   
    prmacro(i)=trSPK(i)/3;
    %prsx(i)=trREL(i)/3;
    
    Temp_inv1=devSPKM(i)^2;
    %Temp_inv2=devREL(i)^2;
    %Stress(i)=sqrt(Temp_inv1+Temp_inv2); 
    stressmacro(i)=devSPKM(i); 

end




 prsmacro=linspace(0,-10e4);
 prx=linspace(0,-10e4);
%  devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
%  devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
%  YieldSurface0standard=Aphi0-Bphi0*pr+Aphi0x-Bphi0x*prx;
%  w=Aphi-Bphi*pr+Aphix-Bphix*prx;
  YieldSurface0macro=abs(Aphi0-Bphi0*prsmacro);
  YieldSurfacemacro=Aphi-Bphi*prsmacro;






figure(1)
plot(-e33standard,-s33standard,'-b',-e33Macro,-s33Macro,'k',-e33comb,-s33comb,'-r')
xlabel('-e33','FontSize',18)
ylabel('-s33','FontSize',18)
 h=legend('Standard plasticity','Macro plasticity','Combined plasticity','Location','SouthEast');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',32)




figure(2)
plot(Bphicomb*prcomb+Bphixcomb*prsxcomb,stresscomb,'-+r',prmacro,stressmacro,'-k',prstandard,stressstandard,'-b',Bphicomb*prscomb+Bphixcomb*prxcomb,YieldSurfacecomb,'-.r',prsmacro,YieldSurfacemacro,'-.k',prsstandard,YieldSurfacestandard,'-.b')
%xlabel('time','FontSize',18)
%ylabel('devrel','FontSize',18)
 h=legend('Combined plasticity','Macro plasticity','Standard plasticity','Fcomb','Fmacro','Fstandard','Location','SouthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',32)




figure(3)
plot(times,devSPKs,'-b',timeM,devSPKM,'k',timecomb,devSPKcomb,'-r')
xlabel('time','FontSize',18)
ylabel('devSPK','FontSize',18)
 h=legend('Standard plasticity','Macro plasticity','Combined plasticity','Location','SouthEast');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',32)

figure(4)
plot(prmacro,devSPKM,'k',prstandard,devSPKs,'-b',prcomb,devSPKcomb,'-r')
xlabel('pressure','FontSize',18)
ylabel('devSPK','FontSize',18)
 h=legend('Macro plasticity','Standard plasticity','Combined plasticity','Location','SouthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',32)


