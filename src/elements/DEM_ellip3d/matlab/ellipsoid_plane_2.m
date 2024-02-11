syms a b c d e f g h i j p q r s x y z real
A=[2*a d f p; d 2*b e q; f e 2*c r; p q r 0];
B=[-g -h -i -s]';
X=A\B;
syms sub;
[Y,sub]=subexpr(X,sub);
ccode(X)
ccode(Y)
ccode(sub)