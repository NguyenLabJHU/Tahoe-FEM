k=1.4; % gas ratio of specific heat capacity
Ma=[1 1.5 1.66 2 5] % 1:0.25:5; [1 1.66 2 5]; freestream Mach number
Ma2=((k-1).*Ma.^2+2)./(2*k.*Ma.^2-(k-1))
P2ToP1 = 2*k/(k+1)*Ma.^2 - (k-1)/(k+1)
P2starToP1= ((k+1)/2*Ma.^2) .^ (k/(k-1)) .* (2*k/(k+1) .* Ma.^2 - (k-1)/(k+1)) .^ (1/(1-k))

D=1; % radius of the blunt obstacle
rho2_to_rho1 = (k+1)*Ma.^2 ./ (2+(k-1)*Ma.^2)
delta=1.1*D./rho2_to_rho1
plot(Ma,P2ToP1,'b-*',Ma,P2starToP1,'r-d', 'LineWidth', 2);
xlabel('free stream Mach number');
ylabel('ratio of downstream pressure to upstream pressure');
legend('pressure behind shock', 'stagnation pressure', 'location', 'best');
