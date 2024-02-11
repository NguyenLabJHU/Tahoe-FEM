% Air
Tc=132.41;
Pc=37.25*1.01325E+5;
R=8.31;
T=1500;
p=1.5E+10; % 15GPa
a=27*R^2*Tc^2/(64*Pc)
b=R*Tc/(8*Pc)
coef=[1 -b-R*T/p a/p -a*b/p];
rt=roots(coef)
moleVol=rt(1);
moleMass=44;
density=(moleMass*1.0E-3)/moleVol % Kg/m^3
