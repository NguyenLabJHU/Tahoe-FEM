config = 3; 
% 1: for spherical wave w/o particles
% 2: for spherical wave impacting thousands of particle from inside a cavity, 3000 particles
% 3: for spherical wave impacting thousands of particle from inside a cavity, 18000 particles
% 3: for spherical wave impacting thousands of particle from inside a cavity, 48000 particles
if config == 1
  z1Coord= 3.5e-2;
  z2Coord= 5.0e-2;
  x1Coord=-1.25e-3;
  x2Coord= 1.25e-3;
elseif config == 2
  z1Coord= 1.3e-1;
  z2Coord= 2.5e-1;
  x1Coord= 0.99e-1;
  x2Coord= 1.01e-1;
elseif config == 3
  z1Coord= 1.0e-1;
  z2Coord= 5.0e-1;
  x1Coord= 1.495e-1;
  x2Coord= 1.51e-1;
elseif config == 4
  z1Coord= 1.0e-1;
  z2Coord= 6.5e-1;
  x1Coord= 1.99e-1;
  x2Coord= 2.01e-1;
end
y1Coord= x1Coord;
y2Coord= x2Coord;

newData1 = importdata('couple_fluidplot_000.dat', ' ', 2); %delimiterIn,headerlinesIn
%disp(newData1);
numeric  = newData1.data;

exNumeric = numeric(( numeric(:, 1) >= x1Coord & numeric(:, 1) <= x2Coord & ...
                      numeric(:, 2) >= y1Coord & numeric(:, 2) <= y2Coord & ...
                      numeric(:, 3) >= z1Coord & numeric(:, 3) <= z2Coord ), :);

fprintf(1, '%15s%15s%15s\n', 'coordX', 'coordY', 'coordZ');
for i = 1:size(exNumeric,1)
  fprintf(1, '%15.6e%15.6e%15.6e\n', exNumeric(i, 1), exNumeric(i, 2), exNumeric(i, 3));
end


