clear all
close all
%
sim=load('eldataloc_1.txt');
%end_step=length(sim(:,1))
end_step=27
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
loc_flag=(sim(1:end_step,24));

I1=sig11+sig22+sig33;

sig1=sig22;
sig2=sig11;
sig3=sig33;
sig1_x = (1/sqrt(2))*(sig1-sig3);
sig1_y = sqrt(2/3)*(sig2-0.5*sig3-0.5*sig1);

xi1=sig22-alpha22;
xi2=sig11-alpha11;
xi3=sig33-alpha33;
xi1_x = (1/sqrt(2))*(xi1-xi3);
xi1_y = sqrt(2/3)*(xi2-0.5*xi3-0.5*xi1);

% A = 843.02;  %#A - failure surface fitting parameters A,B,C
% B = 2.731e-4; % #B
% C = 821.92; %#C

A = 689.2; %  #A - failure surface fitting parameters A,B,C
B = 3.94e-4; % #B
C = 675.2; % #C
N = 6.0;
theta = 0.0;
R=28;    % unitless

%radius = (2/3)*sqrt(2)*(A - C*exp(B*I1(1)) - theta*I1(1));
center = 0;
handle =1;
psi = .8;

figure(1);
%plot(sig1_x, sig1_y, '-o');
plot(xi1_x, xi1_y, '-.o','LineWidth',1);
hold on;

kappa1=-8.05
X=kappa1-R*(A-C*exp(B*kappa1)-theta*kappa1)
I1use=min(I1);
%Fc=1+0.5*(min(I1)-kappa1)*(abs(min(I1)-kappa1)-(min(I1)-kappa1))/((X-kappa1)^2)
Fc=1+0.5*(I1use-kappa1)*(abs(I1use-kappa1)-(I1use-kappa1))/((X-kappa1)^2)
factor=sqrt(2);
radius = factor*sqrt(Fc)*(A - C*exp(B*I1use) - theta*I1use - N)
%center = sig1_x(18) + (sqrt(3)/2)*radius + i*(sig1_y(18) - radius/2)
yield = plot_gamma_axes(center, radius, psi, handle); 

kappa_last=kappa(end_step)
X=kappa_last-R*(A-C*exp(B*kappa_last)-theta*kappa_last)
Fc=1+(I1(end_step)-kappa_last)*(abs(I1(end_step)-kappa_last)-(I1(end_step)-kappa_last))/(2*(X-kappa_last)^2)
radius = factor*sqrt(Fc)*(A - C*exp(B*I1(end_step)) - theta*I1(end_step) - N)
%center = min(sig1_x) + (sqrt(3)/2)*radius + i*(max(sig1_y) - radius/2);
%Gamma_max = 1/psi;
%center = sig1_x(end_step) - (sqrt(3)/2)*radius/Gamma_max + i*(sig1_y(end_step) + radius/(2*Gamma_max))
%center = sig1_x(end_step) - (sqrt(3)/2)*radius/Gamma_max + i*(sig1_y(end_step) + radius/(2*Gamma_max))
final_yield = plot_gamma_axes(center, radius, psi, handle);
