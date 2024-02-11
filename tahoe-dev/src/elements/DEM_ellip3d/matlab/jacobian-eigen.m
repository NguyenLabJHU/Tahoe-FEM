% g: gamma of air
% a: sound speed
% H: total specific enthalpy

% solution of original matrix
clear;clc;
syms r a u v w H g;
r = g - 1;
H = 1/2*(u^2+v^2+w^2)+a^2/(g-1);
A= [    0,                  1,          0,      0,      0;
        r*H-u^2-a^2,        (3-g)*u,    -r*v,   -r*w,   r;
        -u*v,               v,          u,      0,      0;
        -u*w,               w,          0,      u,      0;
        u*((g-2)*H-a^2),    H-r*u^2,    -r*u*v, -r*u*w, g*u];
[VV,D]=eig(A)

v1=VV(:,1)
v2=VV(:,2)
v3=VV(:,3)
v4=VV(:,4)
v5=VV(:,5)
A*v1
u*v1
jordan(A);
% linear combinations are still solutions, but they can simplify formulations greatly
v1+v3*v
v2+v3*w
simplify(v1*v+v2*w+v3*(u^2+v^2+w^2)/2)
simplify(v4/v4(1:1))
simplify(v5/v5(1:1))

% solution of modified matrix with porosity
% f: f = porosity
clear;clc;
syms a u v w H g r f;
r = g - 1;
H = 1/2*(u^2+v^2+w^2)+a^2/(g-1);
A= [    0,     1/f,      0,     0,     0;
      r*f*H-a^2/f-u^2, (2*f-g+1)/f*u, -1/f*r*v, -1/f*r*w,1/f*r;
        -u*v,    v,      u,      0,     0;
        -u*w,   w,     0,      u,     0;
      (r-f)/f*u*(H-a^2/r)-(r+f)/f*u*a^2/g/r, H-r/f*u^2+(1-f)/f/g*a^2, -r/f*u*v, -r/f*u*w, (f+r)/f*u ] ;
[VV,D]=eig(A)

v1=VV(:,1)
v2=VV(:,2)
v3=VV(:,3)
v4=VV(:,4)
v5=VV(:,5)
simplify(v1/v1(1:1))
simplify(v2/v2(1:1))
v3+v5*v
v4+v5*w
simplify(v3*v+v4*w+v5*(u^2+v^2+w^2)/2)
ccode(D(1,1))

