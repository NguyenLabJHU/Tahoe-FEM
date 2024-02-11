syms x y z X Y Z x0 y0 z0 a b c l1 m1 n1 l2 m2 n2 l3 m3 n3 real;
Q = [l1 l2 l3; m1 m2 m3; n1 n2 n3];
localfun = X^2/a^2 + Y^2/b^2 + Z^2/c^2 - 1;
B = [x-x0 y-y0 z-z0]';
vec = Q.'*B;
globalfun = expand(subs(localfun, {X, Y, Z}, {vec(1,1),vec(2,1),vec(3,1)}));

syms H
a1 = collect( subs( subs(globalfun, x^2, H), {x y z}, {0 0 0} ) , H)