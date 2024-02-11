close all;
clear all;
lambda=28.85e6;mu=19.23e6;kf=1e-5; % kf means hydraulic conductivity
ns=0.58;nf=0.42;rhof=1000;rhos=2700;
Sv=nf^2*rhof*9.81/kf;
a=(ns^2*nf*rhof+nf^2*ns*rhos)/(lambda+2*mu)/(nf^2);
b=Sv/(lambda+2*mu)/(nf^2);
traction=4e4;

data1=load('column20m_dynamic_run130_362.txt'); % top
data2=load('column20m_dynamic_run130_344.txt'); % 1m below
data3=load('column20m_dynamic_run130_326.txt'); % 2m below
data4=load('column20m_dynamic_run130_312.txt'); % 3.5m below

time=data1(:,1);

d_y1=data1(:,2);
d_y2=data2(:,2);
d_y3=data3(:,2);
d_y4=data4(:,2);

s22_1=data1(:,3);
s22_2=data2(:,3);
s22_3=data3(:,3);
s22_4=data4(:,3);

p1=data1(:,4);
p2=data2(:,4);
p3=data3(:,4);
p4=data4(:,4);

z=0;
dt=0.0004;
for k=1:1001
  t(k)=(k-1)*dt;
  tau=linspace(0,t(k),1000);
  dtau=tau(2)-tau(1);
  u1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.* ...
                                                    sqrt(tau.^2-a*z^2)./(2*a)).*heaviside(tau-sqrt(a)*z);
  u1(1)=0;
  u(k)=-1/sqrt(a)/(lambda+2*mu)*(sum(u1)*dtau);
  q1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*tau./(2*a));
  q(k)=-1/sqrt(a).*(sum(q1)*dtau);
  g(k)=1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*sqrt(t(k)^2-a*z^2)/(2* ...
                                                   a))*heaviside(t(k)-sqrt(a)*z)-1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*t(k)/(2*a));  
  sigma1=0;
  sigma2(k)=-traction*(1-cos(50.*(t(k)-sqrt(a)*z)))*heaviside(t(k)-sqrt(a)*z).*exp(-b/2/sqrt(a)*z);
  sigma(k)=b/2/sqrt(a)*(sum(sigma1)*dtau)+sigma2(k);
end

figure(1);
plot(t,-u.*1000,'b')
hold on;
plot(time,d_y1.*1000,'bo')
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(mm)')

for k=1:1001
  l(k)=0;
  for k1=1:k
    l(k)=l(k)+q(k-k1+1)*g(k1)*dt;
  end
end

l_dot = diff(l)./diff(t);
l_dot2 = diff(l_dot)./diff(t(1:end-1));
p=1/nf^2/(lambda+2*mu).*(ns*nf*rhof.*l_dot2+Sv.*l_dot(1:end-1));

p=zeros(1,999);
figure(2);
plot(t(1:end-2),-p./1000,'b')
hold on;
plot(time,p1./1000,'bo')
xlabel('TIME(seconds)')
ylabel('EXCESS PORE PRESSURE(kPa)')

figure(3);
plot(t,sigma./1000,'b')
hold on;
plot(time,s22_1./1000,'bo')
xlabel('TIME(seconds)')
ylabel('EFFECTIVE STRESS(kPa)')

z=1.0;
dt=0.0004;
for k=1:1001
  t(k)=(k-1)*dt;
  tau=linspace(0,t(k),1000);
  dtau=tau(2)-tau(1);
  u1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*sqrt(tau.^2-a*z^2)./(2*a)).*heaviside(tau-sqrt(a)*z);
  u(k)=-1/sqrt(a)/(lambda+2*mu)*(sum(u1)*dtau);
  q1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*tau./(2*a));
  q(k)=-1/sqrt(a).*(sum(q1)*dtau);
  g(k)=1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*sqrt(t(k)^2-a*z^2)/(2* ...
                                                   a))*heaviside(t(k)-sqrt(a)*z)-1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*t(k)/(2*a));  
  sigma1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(1,b.* ...
                                                    sqrt(tau.^2-a*z^2)./(2*a)).*z./sqrt(tau.^2-a*z^2).*heaviside(tau-sqrt(a)*z);
  sigma2(k)=-traction*(1-cos(50.*(t(k)-sqrt(a)*z)))*heaviside(t(k)-sqrt(a)*z).*exp(-b/2/sqrt(a)*z);
  sigma(k)=b/2/sqrt(a)*(sum(sigma1)*dtau)+sigma2(k);
end

figure(1);
plot(t,-u.*1000,'r')
hold on;
plot(time,d_y2.*1000,'ro')
legend('y=1.0m')

for k=1:1001
  l(k)=0;
  for k1=1:k
    l(k)=l(k)+q(k-k1+1)*g(k1)*dt;
  end
end

l_dot = diff(l)./diff(t);
l_dot2 = diff(l_dot)./diff(t(1:end-1));
p=1/nf^2/(lambda+2*mu).*(ns*nf*rhof.*l_dot2+Sv.*l_dot(1:end-1));

figure(2);
plot(t(1:end-2),-p./1000,'r')
hold on;
plot(time,p2./1000,'ro')

figure(3);
plot(t,sigma./1000,'r')
hold on;
plot(time,s22_2./1000,'ro')

z=2.0;
dt=0.0004;
for k=1:1001
  t(k)=(k-1)*dt;
  tau=linspace(0,t(k),1000);
  dtau=tau(2)-tau(1);
  u1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*sqrt(tau.^2-a*z^2)./(2*a)).*heaviside(tau-sqrt(a)*z);
  u(k)=-1/sqrt(a)/(lambda+2*mu)*(sum(u1)*dtau);
  q1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*tau./(2*a));
  q(k)=-1/sqrt(a).*(sum(q1)*dtau);
  g(k)=1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*sqrt(t(k)^2-a*z^2)/(2* ...
                                                   a))*heaviside(t(k)-sqrt(a)*z)-1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*t(k)/(2*a));  
  sigma1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(1,b.* ...
                                                    sqrt(tau.^2-a*z^2)./(2*a)).*z./sqrt(tau.^2-a*z^2).*heaviside(tau-sqrt(a)*z);
  sigma2(k)=-traction*(1-cos(50.*(t(k)-sqrt(a)*z)))*heaviside(t(k)-sqrt(a)*z).*exp(-b/2/sqrt(a)*z);
  sigma(k)=b/2/sqrt(a)*(sum(sigma1)*dtau)+sigma2(k);
end

figure(1);
plot(t,-u.*1000,'k')
hold on;
plot(time,d_y3.*1000,'ko')
legend('y=2.0m')

for k=1:1001
  l(k)=0;
  for k1=1:k
    l(k)=l(k)+q(k-k1+1)*g(k1)*dt;
  end
end

l_dot = diff(l)./diff(t);
l_dot2 = diff(l_dot)./diff(t(1:end-1));
p=1/nf^2/(lambda+2*mu).*(ns*nf*rhof.*l_dot2+Sv.*l_dot(1:end-1));

figure(2);
plot(t(1:end-2),-p./1000,'k')
hold on;
plot(time,p3./1000,'ko')

figure(3);
plot(t,sigma./1000,'k')
hold on;
plot(time,s22_3./1000,'ko')

z=3.5;
dt=0.0004;
for k=1:1001
  t(k)=(k-1)*dt;
  tau=linspace(0,t(k),1000);
  dtau=tau(2)-tau(1);
  u1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*sqrt(tau.^2-a*z^2)./(2*a)).*heaviside(tau-sqrt(a)*z);
  u(k)=-1/sqrt(a)/(lambda+2*mu)*(sum(u1)*dtau);
  q1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(0,b.*tau./(2*a));
  q(k)=-1/sqrt(a).*(sum(q1)*dtau);
  g(k)=1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*sqrt(t(k)^2-a*z^2)/(2* ...
                                                   a))*heaviside(t(k)-sqrt(a)*z)-1/sqrt(a)*exp(-b/2/a*t(k))*besseli(0,b*t(k)/(2*a));  
  sigma1=-traction*(1-cos(50.*(t(k)-tau))).*exp(-b/2/a.*tau).*besseli(1,b.* ...
                                                    sqrt(tau.^2-a*z^2)./(2*a)).*z./sqrt(tau.^2-a*z^2).*heaviside(tau-sqrt(a)*z);
  sigma2(k)=-traction*(1-cos(50.*(t(k)-sqrt(a)*z)))*heaviside(t(k)-sqrt(a)*z).*exp(-b/2/sqrt(a)*z);
  sigma(k)=b/2/sqrt(a)*(sum(sigma1)*dtau)+sigma2(k);
end

figure(1);
plot(t,-u.*1000,'m')
hold on;
plot(time,d_y4.*1000,'mo')
legend('ANALYTICAL','FEM')

for k=1:1001
  l(k)=0;
  for k1=1:k
    l(k)=l(k)+q(k-k1+1)*g(k1)*dt;
  end
end

l_dot = diff(l)./diff(t);
l_dot2 = diff(l_dot)./diff(t(1:end-1));
p=1/nf^2/(lambda+2*mu).*(ns*nf*rhof.*l_dot2+Sv.*l_dot(1:end-1));

figure(2);
plot(t(1:end-2),-p./1000,'m')
hold on;
plot(time,p4./1000,'mo')
legend('ANALYTICAL','FEM')

figure(3);
plot(t,sigma./1000,'m')
hold on;
plot(time,s22_4./1000,'mo')
legend('ANALYTICAL','FEM')

figure(1);
text(0.2,0,'z=3.5');
text(0.14,-0.7e-1,'z=2');
text(0.1,-1.4e-1,'z=1');
text(0.08,-2.7e-1,'z=0');

figure(2);
text(0.05,0.3,'z=0');
text(0.06,2,'z=1');
text(0.055,5,'z=2');
text(0.1,4,'z=3.5');

figure(3);
text(0.04,-3.5,'z=0');
text(0.07,-2,'z=1');
text(0.08,-0.7,'z=2');
text(0.07,0.2,'z=3.5');

