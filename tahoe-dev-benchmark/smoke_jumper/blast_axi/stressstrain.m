clear all
%
sim=load('nodal-soil-axi_411.txt');
end_step=length(sim(:,1))
%
time=(sim(1:end_step,1));
dx=(sim(1:end_step,4));
dy=(sim(1:end_step,5));
srr=(sim(1:end_step,6));
szz=(sim(1:end_step,7));
srz=(sim(1:end_step,8));
stt=(sim(1:end_step,9));
alpha11=(sim(1:end_step,10));
alpha22=(sim(1:end_step,11));
alpha33=(sim(1:end_step,12));
alpha23=(sim(1:end_step,13));
alpha13=(sim(1:end_step,14));
alpha12=(sim(1:end_step,15));
kappa=(sim(1:end_step,16));

% Stres Invariants
meanstress=(srr + szz + stt)./3;
I1=meanstress.*3;
I2=(srr.*szz-srz.*srz)+(szz.*stt)+(srr.*stt);
for i=1:end_step
    A=[srr(i) srz(i) 0 ;
       srz(i) szz(i) 0 ;
       0   0   stt(i) ];
    I3(i,1)=det(A);
end

% Principal Stresses
for i=1:end_step
    p = [1 -I1(i) I2(i) -I3(i)];
    r = roots(p);
    sig1(i,1) = max(r);
    sig2(i,1) = r(2);
    sig3(i,1) = min(r);
end

% Principal Stresses
% for i=1:end_step
%     sigM = [srr(i) srz(i) 0 ;
%             srz(i) szz(i) 0 ;
%             0   0   stt(i) ];
%     lambda = eig(sigM);
%     sig1a(i,1) = max(lambda);
%     sig2a(i,1) = lambda(2);
%     sig3a(i,1) = min(lambda);
% end
% 
% I1a=sig1a+sig2a+sig3a;

% Stress Deviator
s1 = sig1-meanstress;
s2 = sig2-meanstress;
s3 = sig3-meanstress;

% J2
sigJ2 = 0.5*(s1.*s1 + s2.*s2 + s3.*s3);



% Plot Failure Surface

A=0.446;  % MPa
B=0; % 1/MPa
C=0;  %MPa
R=10;    % unitless
theta=0.00053; % radians
N=.2;     % MPa
D1=3e-1;
D2=0.0;
W=.078;


kappa0=-500 %MPa
X0=kappa0-R*(A-C*exp(B*kappa0)-theta*kappa0)
kappa_last=(sim(length(time),16))
X=kappa_last-R*(A-C*exp(B*kappa_last)-theta*kappa_last)
I1s=X:.1:100;
%shear failure surface
Ff=A-C*exp(B*I1s)-theta*I1s;
Ff0=(A-N)-C*exp(B*I1s)-theta*I1s;
%Ff=A-C*exp(B*I1);
Ff=0.5*(abs(Ff)+Ff);
Ff0=0.5*(abs(Ff0)+Ff0);
%cap surface
Fc0=1+(I1s-kappa0).*(abs(I1s-kappa0)-(I1s-kappa0))/(2*(X0-kappa0)^2);
Fc0=0.5*(abs(Fc0)+Fc0);
Fc=1+(I1s-kappa_last).*(abs(I1s-kappa_last)-(I1s-kappa_last))/(2*(X-kappa_last)^2);
Fc=0.5*(abs(Fc)+Fc);
%Fy=0.5*(abs(Ff-N)+Ff-N);
Fy=0.5*(abs(Ff)+Ff);
Fy0=0.5*(abs(Ff0)+Ff0);
fJ2yield0=Fc0.*Fy.^2;
fJ2yieldsq0=sqrt(Fc0).*Fy0;
fJ2failsq0=sqrt(Fc0).*Ff;
fJ2yieldsq=sqrt(Fc).*Fy;
fJ2failsq=sqrt(Fc).*Ff;


plot(I1,sqrt(sigJ2),'o',I1s,fJ2yieldsq)