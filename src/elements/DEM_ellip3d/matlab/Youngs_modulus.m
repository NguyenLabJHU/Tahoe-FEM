% reference: A physically based approach to granular media mechanics, David Cole.
% fig 12, crushed gneiss, monotonic loading
P=20
v=0.25
R=1.0e-3
disp=3.25e-5
E=3/sqrt(2)*P*(1-v^2)/(R^0.5*disp^1.5) = 6.7886e+09

% fig 12, weathered aggregate, monotonic loading
P=14
v=0.25
R=1.67e-3
disp=4.25e-5
E=3/sqrt(2)*P*(1-v^2)/(R^0.5*disp^1.5) = 2.4590e+09

% fig 14, crushed gneissgneiss (Nor1), sinusoidal loading
P=12
v=0.25
R=1.0e-3
disp=0.85e-5
E=3/sqrt(2)*P*(1-v^2)/(R^0.5*disp^1.5) = 3.0453e+10

% fig 15, ball-milled gneiss (Nor3), sinusoidal loading
P=14
v=0.25
R=1.41e-3
disp=1.0e-5
E=3/sqrt(2)*P*(1-v^2)/(R^0.5*disp^1.5) = 2.3447e+10
