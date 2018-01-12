rho=1.6e+3;
g=9.8;
h=0:10:200;
stress=rho*g*h/1000;
plot(h,stress,'b-*');
xlabel('depth (meter)');
ylabel('max stress (kPa)');
% 20m  - 320kPa
% 200m - 3.2MPa
% 6cm  - 1kPa
% 50cm - 8kPa
