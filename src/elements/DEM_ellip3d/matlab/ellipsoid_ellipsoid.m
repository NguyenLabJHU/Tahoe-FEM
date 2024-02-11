syms a1 b1 c1 d1 e1 f1 g1 h1 i1 j1 real;
syms a2 b2 c2 d2 e2 f2 g2 h2 i2 j2 real;
syms x y z real;
syms lam real;
f1 = a1*x^2 + b1*y^2 + c1*z^2 + d1*x*y + e1*y*z + f1*z*x + g1*x + h1*y + i1*z + j1;
f2 = a2*x^2 + b2*y^2 + c2*z^2 + d2*x*y + e2*y*z + f2*z*x + g2*x + h2*y + i2*z + j2;
F = f1 + lam*f2;
eq1 = diff(F,x);
eq2 = diff(F,y);
eq3 = diff(F,z);
eq4 = diff(F,lam);
S = solve(eq1,eq2,eq3,x,y,z); % system of linear equations
eq4 = expand(subs(eq4,{x y z},{S.x S.y S.z}));
out = collect(eq4,lam);
fid = fopen('ellipsoid_ellipsoid.txt','w');
fprintf(fid,'%s\n',char(out));
fclose(fid);

