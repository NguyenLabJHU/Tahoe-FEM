sigma=0:1.0e+2:1.0e+4;
E=2.5e+10;
ratio=1.875*(sigma ./ E) .^(2.0/3);
plot(sigma, ratio, 'b-*');
xlabel('stress (Pa)');
ylabel('static penetr ratio');
