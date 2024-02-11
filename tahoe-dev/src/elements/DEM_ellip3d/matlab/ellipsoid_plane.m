% plane
syms x y z real;
syms p q r s real;
plane =  p*x + q*y + r*z + s;

% ellipsoid
syms a b c d e f g h i j real;
syms lam real;
ellipsoid = a*x^2 + b*y^2 + c*z^2 + d*x*y + e*y*z + f*z*x + g*x + h*y + i*z + j;
F = ellipsoid + lam*plane;
eq1=diff(F,x);
eq2=diff(F,y);
eq3=diff(F,z);
eq4=diff(F,lam);
S = solve(eq1,eq2,eq3,eq4,x,y,z,lam); % system of linear equations

% print
fid = fopen('ellipsoid_plane.txt','w');
fprintf(fid,'%s\n',char(ccode(S.x)));
fprintf(fid,'%s\n',char(ccode(S.y)));
fprintf(fid,'%s\n',char(ccode(S.z)));
fclose(fid);
