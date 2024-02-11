clear all
%

sim=load('eldata_restart_1.txt');
test5 = importdata('sj_test_5.csv',',',9);
test5_data = test5.data;
test6 = importdata('sj_test_6.csv',',',9);
test6_data = test6.data;
end_step=length(sim(:,1))
%
time=(sim(1:end_step,1));
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
%ip8
alph11=(sim(1:end_step,14));
alph22=(sim(1:end_step,15));
alph33=(sim(1:end_step,16));
alph23=(sim(1:end_step,17));
alph13=(sim(1:end_step,18));
alph12=(sim(1:end_step,19));
kappa=(sim(1:end_step,20));
%press=sim(1:end_step,21);
J2=(sim(1:end_step,22));
J3=(sim(1:end_step,23));
loc_flag=(sim(1:end_step,24));
%
sim=load('eldata_1.txt');
sig11f=(sim(1:end_step,2));
sig22f=(sim(1:end_step,3));
eps22f=(sim(1:end_step,9));
%
% sim=load('eldata_unload_1.txt');
% sig11u=(sim(1:end_step,2));
% sig22u=(sim(1:end_step,3));
% eps22u=(sim(1:end_step,9));
%
press=(sig11+sig22+sig33)/3;
I1s=3*press;
devsig11=sig11-press;
devsig22=sig22-press;
devsig33=sig33-press;
xi11=devsig11-alph11;
xi22=devsig22-alph22;
xi33=devsig33-alph33;
xi23=sig23-alph23;
xi12=sig12-alph12;
xi13=sig13-alph13;
xiprod=xi11.^2+xi22.^2+xi33.^2+2*xi12.^2+2*xi23.^2+2*xi13.^2;
J2xi=.5*xiprod;
sJ2xi=sqrt(J2xi);
sJ2=sqrt(J2);
%
A=0.66;  % MPa
B=0; % 1/MPa
C=0;  %MPa
R=10;    % unitless
theta=0.0003; % radians
N=.12;     % MPa
D1=3e-1;
D2=0.0;
W=.078;

%X0=-460; % MPas
%kappa0=-1.32; %MPa
kappa0=-190 %MPa
X0=kappa0-R*(A-C*exp(B*kappa0)-theta*kappa0)
kappa_last=(sim(length(time),20))
X=kappa_last-R*(A-C*exp(B*kappa_last)-theta*kappa_last)
%I1=-1725:1:50;
I1=X:.1:100;
%shear failure surface
Ff=A-C*exp(B*I1)-theta*I1;
Ff0=(A-N)-C*exp(B*I1)-theta*I1;
%Ff=A-C*exp(B*I1);
Ff=0.5*(abs(Ff)+Ff);
Ff0=0.5*(abs(Ff0)+Ff0);
%cap surface
Fc0=1+(I1-kappa0).*(abs(I1-kappa0)-(I1-kappa0))/(2*(X0-kappa0)^2);
Fc0=0.5*(abs(Fc0)+Fc0);
Fc=1+(I1-kappa_last).*(abs(I1-kappa_last)-(I1-kappa_last))/(2*(X-kappa_last)^2);
Fc=0.5*(abs(Fc)+Fc);
%Fy=0.5*(abs(Ff-N)+Ff-N);
Fy=0.5*(abs(Ff)+Ff);
Fy0=0.5*(abs(Ff0)+Ff0);
fJ2yield0=Fc0.*Fy.^2;
fJ2yieldsq0=sqrt(Fc0).*Fy0;
fJ2failsq0=sqrt(Fc0).*Ff;
fJ2yieldsq=sqrt(Fc).*Fy;
fJ2failsq=sqrt(Fc).*Ff;
%
%
figure(5)
%plot(abs(eps22),abs(sig22)-abs(sig11),'LineWidth',1)
plot(abs(eps22),abs(sig22)-abs(sig11),'x')
hold on
plot(abs(eps22f),abs(sig22f)-abs(sig11f),'x')
%plot(abs(eps22u),abs(sig22u)-abs(sig11u),'+')
plot((test5_data(:,1)/100)+3e-3,test5_data(:,2),'r','LineWidth',1.5);
plot((test6_data(:,1)/100)+3e-3,test6_data(:,2),'b','LineWidth',1.5);
xlabel('STRAIN')
ylabel('STRESS (MPa)')
%legend('100 steps','40 steps')
set(gca,'FontName','Helvetica','FontSize',16)
%
%
figure(9)
%plot(I1s,sJ2xi,'o',I1,fJ2yieldsq0,I1,fJ2yieldsq2,I1,fJ2yieldsq1,I1,fJ2yieldsq,'LineWidth',1)
plot(I1s,sJ2,'o',I1,fJ2yieldsq0,I1,fJ2yieldsq,'LineWidth',1)
xlabel('I1')
ylabel('sJ2xi')
legend(...
'stress path',...
'initial yield surface',...
'final yield surface')
%'intermediate yield surface',...
%'intermediate yield surface',...
%'final yield surface')
set(gca,'FontName','Helvetica','FontSize',16)
%