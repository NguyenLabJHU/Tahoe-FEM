% cylinder
syms x y z real;
syms x0 y0 z0 l3 m3 n3 R real; % (x0 y0 z0)-point, (l3 m3 n3)-z direction, R-radius
cylinder = (1-l3)*(x-x0)^2 + (1-m3)*(y-y0)^2 + (1-n3)*(z-z0)^2 - R^2;
cylinder = subs(cylinder,{l3 m3 n3},{0 0 1});

% ellipsoid
syms a b c d e f g h i j real;
syms lam real;
ellipsoid = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*y*z + f*z*x + g*x + h*y + i*z + j;
F = ellipsoid + lam*cylinder;
eq1=diff(F,x);
eq2=diff(F,y);
eq3=diff(F,z);
eq4=diff(F,lam);
S = solve(eq1,eq2,eq3,x,y,z); % system of linear equations
eq4 = expand(subs(eq4,{x y z},{S.x S.y S.z}));
out = collect(eq4,lam);

% print
fid = fopen('ellipsoid_cylinder.txt','w');
fprintf(fid,'%s',char(out));
fclose(fid);