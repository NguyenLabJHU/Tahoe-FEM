
close all
clear all

format short e

%sim=load('el_data_1.txt');
sim=load('el_DP_eta_1.txt');

time=sim(:,1);
%ip1
s11=sim(:,2);
s22=sim(:,3);
s33=sim(:,4);
s12=sim(:,5);
s13=sim(:,6);
s21=sim(:,7);
s23=sim(:,8);
s31=sim(:,9);
s32=sim(:,10);
e11=sim(:,11);
e22=sim(:,12);
e33=sim(:,13);
e12=sim(:,14);
e13=sim(:,15);
e21=sim(:,16);
e23=sim(:,17);
e31=sim(:,18);
e32=sim(:,19);

%final
c=(sim(length(time),20));
Zc=(sim(length(time),21));
hc=(sim(length(time),22));
Delgamma=(sim(length(time),23));
%initial
c0=(sim(1,20));
Zc0=(sim(1,21));
hc0=(sim(1,22));
Delgamma0=(sim(1,23));
%c0=2e5;
mid_step=round(length(time)/2);
cmid=(sim(mid_step,20));
%

phi0=30*pi/180;
psi0=30*pi/180;

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
    %devSinv(i)=norm(devs);
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

pr=linspace(0,-10e5);
devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));

figure(1)
plot(pr,devs0inv,'--k',pr,devsinv,'-k',p,devSinv,'-or');
legend('F0','F','stress path')
xlabel('mean stress')
ylabel('||devS||')
set(gca,'FontName','Helvetica','FontSize',16)




figure(2)
plot(-e33,-s33*1e-6,'-or','linewidth',1);
%legend('F0','F','stress path')
xlabel('e_{33}')
ylabel('\sigma_{33}   (MPa)')
set(gca,'FontName','Helvetica','FontSize',16)
%h=legend('FSE-DP plasticity','FSE micromorphic DP plasticity w \eta and \phi_{33}');
% xlabel('Eulerian strain')
% ylabel('Cauchy stress')

% hold on 
% 
% 
% 
% sim=load('eldata_1.txt');
% time=sim(:,1);
% %ip1
% s11=sim(:,2);
% s22=sim(:,3);
% s33=sim(:,4);
% s12=sim(:,5);
% s13=sim(:,6);
% s21=sim(:,7);
% s23=sim(:,8);
% s31=sim(:,9);
% s32=sim(:,10);
% e11=sim(:,11);
% e22=sim(:,12);
% e33=sim(:,13);
% e12=sim(:,14);
% e13=sim(:,15);
% e21=sim(:,16);
% e23=sim(:,17);
% e31=sim(:,18);
% e32=sim(:,19);
% 
% %final
% c=(sim(length(time),20));
% Zc=(sim(length(time),21));
% hc=(sim(length(time),22));
% Delgamma=(sim(length(time),23));
% %initial
% c0=(sim(1,20));
% Zc0=(sim(1,21));
% hc0=(sim(1,22));
% Delgamma0=(sim(1,23));
% %c0=2e5;
% mid_step=round(length(time)/2);
% cmid=(sim(mid_step,20));
% %
% 
% phi0=30*pi/180;
% psi0=30*pi/180;
% 
% beta=-1; %TCs
% 
% Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0))
% Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0))
% Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0))
% Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0))
% 
% Im=[1 0 0; 0 1 0; 0 0 1];
% % Principal Stresses
% for i=1:length(time)
%     sig = [s11(i) s12(i) s13(i) ;
%            s21(i) s22(i) s23(i) ;
%            s31(i) s32(i) s33(i) ];
%     lambda = eig(sig);
%     s1(i) = max(lambda);
%     s2(i) = lambda(2);
%     s3(i) = min(lambda);
%     p(i)=(sig(1,1)+sig(2,2)+sig(3,3))/3;
%     devs=sig-p(i)*Im;
%     Temp_inv=0;
%     for j=1:3
%         for k=1:3
%             Temp_inv=Temp_inv+devs(j,k)*devs(j,k);
%         end
%     end 
%     %devSinv(i)=norm(devs);
%     devSinv(i)=sqrt(Temp_inv);
% end
% 
% I1=s1+s2+s3;
% meanstress=I1./3;
% 
% % Principal Strains
% for i=1:length(time)
%     eps = [e11(i) e12(i) e13(i) ;
%            e21(i) e22(i) e23(i) ;
%            e31(i) e32(i) e33(i) ];
%     lambda = eig(eps);
%     e1(i) = max(lambda);
%     e2(i) = lambda(2);
%     e3(i) = min(lambda);
% end
% 
% pr=linspace(0,-10e5);
% devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
% devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
% 
% figure(1)
% plot(pr,devs0inv,'--k',pr,devsinv,'-k',p,devSinv,'-r');
% legend('F0','F','stress path')
% xlabel('mean stress')
% ylabel('||devS||')
% set(gca,'FontName','Helvetica','FontSize',16)
% 
% plot(-e33,-s33*1e-6,'-k','lineWidth',2);
% legend('F0','F','stress path')
% xlabel('e33')
% ylabel('s33')
% set(gca,'FontName','Helvetica','FontSize',16)
% 
% hold on
% 
% 
% sim=load('eldata_eta_kappa_1.txt');
% time=sim(:,1);
% %ip1
% s11=sim(:,2);
% s22=sim(:,3);
% s33=sim(:,4);
% s12=sim(:,5);
% s13=sim(:,6);
% s21=sim(:,7);
% s23=sim(:,8);
% s31=sim(:,9);
% s32=sim(:,10);
% e11=sim(:,11);
% e22=sim(:,12);
% e33=sim(:,13);
% e12=sim(:,14);
% e13=sim(:,15);
% e21=sim(:,16);
% e23=sim(:,17);
% e31=sim(:,18);
% e32=sim(:,19);
% 
% %final
% c=(sim(length(time),20));
% Zc=(sim(length(time),21));
% hc=(sim(length(time),22));
% Delgamma=(sim(length(time),23));
% %initial
% c0=(sim(1,20));
% Zc0=(sim(1,21));
% hc0=(sim(1,22));
% Delgamma0=(sim(1,23));
% %c0=2e5;
% mid_step=round(length(time)/2);
% cmid=(sim(mid_step,20));
% %
% 
% phi0=30*pi/180;
% psi0=30*pi/180;
% 
% beta=-1; %TCs
% 
% Aphi=2*sqrt(6)*c*cos(phi0)/(3+beta*sin(phi0))
% Aphimid=2*sqrt(6)*cmid*cos(phi0)/(3+beta*sin(phi0))
% Aphi0=2*sqrt(6)*c0*cos(phi0)/(3+beta*sin(phi0))
% Bphi0=2*sqrt(6)*sin(phi0)/(3+beta*sin(phi0))
% 
% Im=[1 0 0; 0 1 0; 0 0 1];
% % Principal Stresses
% for i=1:length(time)
%     sig = [s11(i) s12(i) s13(i) ;
%            s21(i) s22(i) s23(i) ;
%            s31(i) s32(i) s33(i) ];
%     lambda = eig(sig);
%     s1(i) = max(lambda);
%     s2(i) = lambda(2);
%     s3(i) = min(lambda);
%     p(i)=(sig(1,1)+sig(2,2)+sig(3,3))/3;
%     devs=sig-p(i)*Im;
%     Temp_inv=0;
%     for j=1:3
%         for k=1:3
%             Temp_inv=Temp_inv+devs(j,k)*devs(j,k);
%         end
%     end 
%     %devSinv(i)=norm(devs);
%     devSinv(i)=sqrt(Temp_inv);
% end
% 
% I1=s1+s2+s3;
% meanstress=I1./3;
% 
% % Principal Strains
% for i=1:length(time)
%     eps = [e11(i) e12(i) e13(i) ;
%            e21(i) e22(i) e23(i) ;
%            e31(i) e32(i) e33(i) ];
%     lambda = eig(eps);
%     e1(i) = max(lambda);
%     e2(i) = lambda(2);
%     e3(i) = min(lambda);
% end
% 
% pr=linspace(0,-10e5);
% devs0inv=0.5*((Aphi0-Bphi0*pr)+abs(Aphi0-Bphi0*pr));
% devsinv=0.5*((Aphi-Bphi0*pr)+abs(Aphi-Bphi0*pr));
% 
% % figure(1)
% % plot(pr,devs0inv,'--k',pr,devsinv,'-k',p,devSinv,'-r');
% % legend('F0','F','stress path')
% % xlabel('mean stress')
% % ylabel('||devS||')
% % set(gca,'FontName','Helvetica','FontSize',16)
% % 
% plot(-e33,-s33*1e-6,'-+k','lineWidth',2);
% h=legend('FSE-DP plasticity','FSE micromorphic DP plasticity with \eta and \phi_{33}','FSE micromorphic DP plasticity with \eta,\kappa and \phi_{33}');
% %legend('F0','F','stress path')
% % xlabel('e_{33}')
% % ylabel('s33')
% set(gca,'FontName','Helvetica','FontSize',16)

