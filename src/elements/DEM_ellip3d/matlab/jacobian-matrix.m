clear;clc;
syms rho u v w E;
syms u1 u2 u3 u4 u5;
syms f1 f2 f3 f4 f5;
syms g phi; % g=gama
syms J;

f1=u2/phi;
f2=u2^2/u1+1/phi*(g-1)*(u5-1/2/u1*(u2^2+u3^2+u4^2));
f3=u2*u3/u1;
f4=u2*u4/u1;
f5=u2*u5/u1+1/phi*u2/u1*(g-1)*(u5-1/2/u1*(u2^2+u3^2+u4^2));

J=jacobian([f1; f2; f3; f4; f5], [u1 u2 u3 u4 u5])
J1=simplify( subs(J, {u1,u2,u3,u4,u5}, [rho, rho*u, rho*v, rho*w, E]) )
pretty(J1)

