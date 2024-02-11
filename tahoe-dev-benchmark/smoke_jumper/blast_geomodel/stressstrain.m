clear all
close all
%
% A
sim=load('eldata_155.txt');
%
% B
%sim=load('eldata_126.txt');
%
% C
%sim=load('eldata_002.txt');
%
end_step=length(sim(:,1))
%end_step=200
%
time=(sim(1:end_step,1));
%ip1
sig11=(sim(1:end_step,2));
sig22=(sim(1:end_step,3));
sig33=(sim(1:end_step,4));
sig23=(sim(1:end_step,5));
sig13=(sim(1:end_step,6));
sig12=(sim(1:end_step,7));
eps11=(sim(1:end_step,8));
eps22=(sim(1:end_step,9));
eps33=(sim(1:end_step,10));
eps23=(sim(1:end_step,11));
eps13=(sim(1:end_step,12));
eps12=(sim(1:end_step,13));
alpha11=(sim(1:end_step,14));
alpha22=(sim(1:end_step,15));
alpha33=(sim(1:end_step,16));
alpha23=(sim(1:end_step,17));
alpha13=(sim(1:end_step,18));
alpha12=(sim(1:end_step,19));
kappa=(sim(1:end_step,20));
press=sim(1:end_step,21);
J2=(sim(1:end_step,22));
J3=(sim(1:end_step,23));
%loc_flag=(sim(1:end_step,24));

%press=(sig11+sig22+sig33)/3;

% for i=1:end_step
%     sig23(i,1)=0;
%     sig13(i,1)=0;
% end

% Principal Stresses
for i=1:end_step
    sig = [sig11(i) sig12(i) sig13(i) ;
           sig12(i) sig22(i) sig23(i) ;
           sig13(i) sig23(i) sig33(i) ];
    lambda = eig(sig);
    sig1(i,1) = max(lambda);
    sig2(i,1) = lambda(2);
    sig3(i,1) = min(lambda);
end

I1=sig1+sig2+sig3;
meanstress=I1./3;

% % Principal Stresses
% for i=1:end_step
%     p = [1 -I1(i) I2(i) -I3(i)];
%     r = roots(p);
%     sig1(i,1) = max(r);
%     sig2(i,1) = r(2);
%     sig3(i,1) = min(r);
% end

% Principal BackStresses
for i=1:end_step
    alpha = [alpha11(i) alpha12(i) alpha13(i) ;
             alpha12(i) alpha22(i) alpha23(i) ;
             alpha13(i) alpha23(i) alpha33(i) ];
    lambda = eig(alpha);
    a1(i,1) = max(lambda);
    a2(i,1) = lambda(2);
    a3(i,1) = min(lambda);
end

%I1a=sig1a+sig2a+sig3a;

% Stress Deviator
s1 = sig1-meanstress;
s2 = sig2-meanstress;
s3 = sig3-meanstress;

% Back Stress
devsig11=sig11-meanstress;
devsig22=sig22-meanstress;
devsig33=sig33-meanstress;
xi11=devsig11-alpha11;
xi22=devsig22-alpha22;
xi33=devsig33-alpha33;
xi23=sig23-alpha23;
xi12=sig12-alpha12;
xi13=sig13-alpha13;
xiprod=xi11.^2+xi22.^2+xi33.^2+2*xi12.^2+2*xi23.^2+2*xi13.^2;
J2xi=.5*xiprod;
sJ2xi=sqrt(J2xi);

% J2
sigJ2 = 0.5*(s1.*s1 + s2.*s2 + s3.*s3);
sJ2 = sqrt(sigJ2);

% Offset, evolution of alpha
J2a=0.5.*(a1.^2 + a2.^2 + a3.^2);
sJ2a=sqrt(J2a);

% parameters
A=5.35;  % MPa
B=0; % 1/MPa
C=0;  %MPa
R=10;    % unitless
theta=1.8e-3; % radians
N=4;     % MPa
D1=3e-1;
D2=0.0;
W=.078;

kappa0=-10; %MPa
X0=kappa0-R*(A-C*exp(B*kappa0)-theta*kappa0);
kappa_last=sim(length(time),20);
X=kappa_last-R*(A-C*exp(B*kappa_last)-theta*kappa_last);

I1s=-114:.1:50;

%shear failure surface
Ff=A-C*exp(B*I1s)-theta*I1s;
Ff0=(A-N)-C*exp(B*I1s)-theta*I1s;
Ff=0.5*(abs(Ff)+Ff);
Ff0=0.5*(abs(Ff0)+Ff0);

%cap surface
%initial
Fc0=1+(I1s-kappa0).*(abs(I1s-kappa0)-(I1s-kappa0))/(2*(X0-kappa0)^2);
Fc0=0.5*(abs(Fc0)+Fc0);
%final
Fc=1+(I1s-kappa_last).*(abs(I1s-kappa_last)-(I1s-kappa_last))/(2*(X-kappa_last)^2);
Fc=0.5*(abs(Fc)+Fc);

Fy=0.5*(abs(Ff)+Ff);
Fy0=0.5*(abs(Ff0)+Ff0);
fJ2yield0=Fc0.*Fy.^2;
fJ2yieldsq0=sqrt(Fc0).*Fy0;
fJ2failsq0=sqrt(Fc0).*Ff;
fJ2yieldsq=sqrt(Fc).*Fy;
fJ2failsq=sqrt(Fc).*Ff;

Fyxi=0.5*(abs(Ff-N)+Ff-N);
fJ2xiyieldsq0=sqrt(Fc0).*Fyxi;
fJ2xiyieldsq=sqrt(Fc).*Fyxi;

I1s2=-150:.1:50;
%shear failure surface
Ff=A-C*exp(B*I1s2)-theta*I1s2;
%Ff=A-C*exp(B*I1);
Ff=0.5*(abs(Ff)+Ff);


figure(1)
hold on
plot(I1s,fJ2yieldsq0,'LineWidth',1)
plot(I1s2,Ff,'Color',[1 0 0],'LineWidth',1)
plot(I1,sJ2,'o','Color',[0 0.498 0],'LineWidth',1)
plot(I1s,fJ2yieldsq0+sJ2a(end_step),I1s,sJ2a(end_step)-fJ2yieldsq0,'b','LineWidth',1,'LineStyle','-.')
ylim([0 6])
xlabel('I1 (MPa)','FontName','Times','FontSize',12);
ylabel('sqrt(J_2) (MPa)','FontName','Times','FontSize',12);
set(gca,'FontName','Times','XGrid','on','YGrid','on')
legend1=legend('Initial Yield Surface','Limit Surface','Simulation J2','Final Yield Surface');
set(legend1,'Position',[0.1969 0.4655 0.2777 0.1577],'FontName','Times');


figure(2)
hold on
plot(I1s,fJ2xiyieldsq0,'LineWidth',1)
plot(I1,sJ2xi,'x','Color',[0 0.498 0],'LineWidth',1)
plot(I1s,fJ2xiyieldsq,'b','LineWidth',1,'LineStyle','-.')
ylim([0 2])
xlabel('I1 (MPa)','FontName','Times','FontSize',12);
ylabel('sqrt(J_2^\xi) (MPa)','FontName','Times','FontSize',12);
set(gca,'FontName','Times','XGrid','on','YGrid','on')
legend1=legend('Initial Yield Surface','Simulation J2xi','Final Yield Surface');
set(legend1,'Position',[0.1969 0.4655 0.2777 0.1577],'FontName','Times');


% figure(2)
% plot(time,eps22,'LineWidth',1)
% 
% xlabel('Time (Sec)','FontName','Times');
% ylabel('Axial Strain','FontName','Times');
% set(gca,'FontName','Times','XGrid','on','YGrid','on')
% 
% 
% figure(3)
% plot(time,meanstress,'LineWidth',1)
% 
% xlabel('Time (Sec)','FontName','Times');
% ylabel('Mean Stress (MPa)','FontName','Times');
% set(gca,'FontName','Times','XGrid','on','YGrid','on')

% hold on
% plot(I1s,fJ2yieldsq0,I1s,Ff)
% 
% for i=1:end_step
%     hold on
%     if (loc_flag(i)==1)
%         plot(I1(i),sJ2(i),'ro')
%     else
%         plot(I1(i),sJ2(i),'kx')
%     end
% end

% fig4=figure
% set(gcf,'Position', [200 200 800 400]);
% 
% for i=1:end_step
% 
% 
% hold on
% 
% subplot(1,2,1)
% hold on
% 
% plot(time(i),eps22(i),'o','LineWidth',2,'Color',[0 0.498 0])
% 
% xlim([0 0.004])
% ylim([-0.05 0])
% xlabel('Time (Sec)','FontName','Times','FontSize',12);
% ylabel('Axial Strain','FontName','Times','FontSize',12);
% set(gca,'FontName','Times','XGrid','on','YGrid','on')
% 
% 
% subplot(1,2,2)
% hold on
% plot(I1s,fJ2yieldsq0,'LineWidth',2)
% plot(I1s,Ff,'Color',[1 0 0],'LineWidth',2)
% 
% 
% xlim([-150 50])
% ylim([0 6])
% xlabel('I1 (MPa)','FontName','Times','FontSize',12);
% %ylabel('sqrt(J_2^\xi) (MPa)','FontName','Times','FontSize',12);
% ylabel('sqrt(J_2) (MPa)','FontName','Times','FontSize',12);
% set(gca,'FontName','Times','XGrid','on','YGrid','on')
% 
% 
% 
% plot(I1(i),sJ2(i),'o','Color',[0 0.498 0],'LineWidth',2)
% M(i)=getframe;
% 
% 
% 
% end
