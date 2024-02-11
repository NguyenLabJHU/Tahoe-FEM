config = 2; 
%1: 3000 particles; 
%2: 18000 particles
%3: 50000 particles
if config == 1
  z1Coord= 1.0e-1;
  z2Coord= 2.0e-1;
  x1Coord= 9.0e-2;
  x2Coord= 1.1e-1;
elseif config == 2
  z1Coord= 3.2e-1;
  z2Coord= 3.3e-1;
  x1Coord= 1.3e-1;
  x2Coord= 1.7e-1;

% for vertical direction
%{
  z1Coord= 1.0e-1;
  z2Coord= 3.3e-1;
  x1Coord= 1.43e-1;
  x2Coord= 1.57e-1;

% for horizontal direction, combination of z and x
  z1Coord= 3.2e-1;
  z2Coord= 3.3e-1;

  z1Coord= 2.7e-1;
  z2Coord= 2.8e-1;

  z1Coord= 2.45e-1;
  z2Coord= 2.55e-1;

  z1Coord= 2.2e-1;
  z2Coord= 2.3e-1;

  z1Coord= 1.75e-1;
  z2Coord= 1.85e-1;

  x1Coord= 1.95e-1;
  x2Coord= 2.05e-1;

  x1Coord= 2.45e-1;
  x2Coord= 2.55e-1;
%}
elseif config == 3
% for vertical direction
  z1Coord= 1.0e-1;
  z2Coord= 4.4e-1;
  x1Coord= 1.94e-1;
  x2Coord= 2.06e-1;

%{
% for horizontal direction, combination of z and x
  z1Coord= 4.35e-1;
  z2Coord= 4.4e-1;

  z1Coord= 3.12e-1;
  z2Coord= 3.18e-1;

  z1Coord= 2.1e-1;
  z2Coord= 2.2e-1;

  z1Coord= 1.2e-1;
  z2Coord= 1.3e-1;

  x1Coord= 2.5e-1;
  x2Coord= 2.7e-1;

  x1Coord= 3.0e-1;
  x2Coord= 3.2e-1;

  x1Coord= 3.45e-1;
  x2Coord= 3.6e-1;
%}
end

y1Coord= x1Coord;
y2Coord= x2Coord;

newData1 = importdata('ini_particle', ' ', 2); %delimiterIn,headerlinesIn
%disp(newData1);
numeric  = newData1.data;

exNumeric = numeric(( numeric(:, 6) >= x1Coord & numeric(:, 6) <= x2Coord & ...
                      numeric(:, 7) >= y1Coord & numeric(:, 7) <= y2Coord & ...
                      numeric(:, 8) >= z1Coord & numeric(:, 8) <= z2Coord ), :);

fprintf(1, '%10s%15s%15s%15s\n', 'particleID', 'coordX', 'coordY', 'coordZ');
for i = 1:size(exNumeric,1)
  fprintf(1, '%10i%15.6e%15.6e%15.6e\n', exNumeric(i, 1), exNumeric(i, 6), exNumeric(i, 7), exNumeric(i, 8));
end


