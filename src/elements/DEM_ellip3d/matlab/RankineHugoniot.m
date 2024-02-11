rhoR=1.225;
pR=1.01325e+5;
uR=0;
MachShock= 1:1:20

gamma=1.4;
Rs=287;
soundSpeedR=sqrt(gamma*pR/rhoR);
MachR=uR/soundSpeedR;
TempR=soundSpeedR^2/(gamma*Rs)-273.15;
%TempR=pR/(rhoR*Rs)-273.15;

shockSpeed = MachShock * sqrt(gamma*pR/rhoR);
rhoL = (rhoR*(shockSpeed-uR)).^2 * (1+gamma) ./ ( rhoR*(shockSpeed-uR).^2 * (gamma-1) + 2*pR*gamma);
pL   = (pR*(1-gamma)+2*rhoR*(shockSpeed-uR) .^ 2) / (1+gamma)
uL   = ( rhoR*(shockSpeed-uR) .* (2*shockSpeed + uR*(gamma-1)) - 2*pR*gamma ) ./ (rhoR * (shockSpeed-uR) * (1+gamma))

soundSpeedL=sqrt(gamma*pL ./ rhoL);
MachL=uL ./ soundSpeedL
TempL=soundSpeedL .^2/(gamma*Rs)-273.15;
plot(MachShock, MachL);
xlabel('Mach of shock');
ylabel('MachL');
