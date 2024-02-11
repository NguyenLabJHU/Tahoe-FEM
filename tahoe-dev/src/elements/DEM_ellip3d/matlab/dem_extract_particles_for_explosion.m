function dem_extract_particles(fileToRead1)

x1Coord= 3.5e-3;
x2Coord= 4.5e-3;
y1Coord= 9.5e-3;
y2Coord= 10e-3;
z1Coord= 12e-3;
z2Coord= 13e-3;

newData1 = importdata(fileToRead1, ' ', 2); %delimiterIn,headerlinesIn
%disp(newData1);
numeric  = newData1.data;

exNumeric = numeric(( numeric(:, 6) >= x1Coord & numeric(:, 6) <= x2Coord & ...
                      numeric(:, 7) >= y1Coord & numeric(:, 7) <= y2Coord & ...
                      numeric(:, 8) >= z1Coord & numeric(:, 8) <= z2Coord ), :);

fprintf(1, '%10s%15s%15s%15s\n', 'particleID', 'coordX', 'coordY', 'coordZ');
for i = 1:size(exNumeric,1)
  fprintf(1, '%10i%15.6e%15.6e%15.6e\n', exNumeric(i, 1), exNumeric(i, 6), exNumeric(i, 7), exNumeric(i, 8));
end


