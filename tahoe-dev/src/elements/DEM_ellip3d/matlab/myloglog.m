x=0:100; 
y=log(x);
%y=x;
figure
plot(x,y)
figure
semilogx(x,y)
figure
loglog(x,y)

x=0.1:0.1:100
figure
y=log10(x);
plot(log10(x),log10(y))

