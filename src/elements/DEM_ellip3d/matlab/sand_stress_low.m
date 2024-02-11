rho=1.6e+3;
g=9.8;
h=0:0.05:1;
stress=rho*g*h/1000;
plot(h,stress,'b-*');
xlabel('depth (meter)');
ylabel('max stress (kPa)');

