% CO2
Tc=304.19;
Pc=7.38E+6;
R=8.31;
T=1500;
p=1.5E+10; % 15GPa
a_=27*R^2*Tc^2/(64*Pc)
b_=R*Tc/(8*Pc)
a=364.77 * 1.0E+3 * (1.0E-3)^2
b=0.043 * 1.0E-3
coef=[1 -b-R*T/p a/p -a*b/p];
rt=roots(coef)
moleVol=rt(1);
moleMass=44;
density=(moleMass*1.0E-3)/moleVol % Kg/m^3
