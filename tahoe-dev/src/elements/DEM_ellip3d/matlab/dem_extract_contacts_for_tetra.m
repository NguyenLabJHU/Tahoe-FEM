function dem_extract_contacts_for tetra(fileToRead1)

x1Coord= 2.48e-1;
x2Coord= 2.52e-1;
%y1Coord= 9.5e-3;
%y2Coord= 10e-3;
%z1Coord= 12e-3;
%z2Coord= 13e-3;

fIn = fopen(fileToRead1);
text = fscanf(fIn, '%s', 20); % 1 + 16 + 3 = 20 strings in the data file
numeric = fscanf(fIn, '%e', [16 inf]);
fclose(fIn);
numeric = numeric.';
exNumeric = numeric(( numeric(:, 1) >= x1Coord & numeric(:, 1) <= x2Coord ), :);

fOut=fopen(strcat('ex_', fileToRead1), 'w');
fprintf(fOut, '%s\n%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n', ...
              'VARIABLES=', 'x', 'y', 'z', 'normal_x', 'normal_y', 'normal_z', 'tangt_x', 'tangt_y', 'tangt_z', ...
              'total_x', 'total_y', 'total_z', 'normal', 'shear', 'total', 'penetration');
fprintf(fOut, 'ZONE I=%d, DATAPACKING=POINT\n', size(exNumeric,1));
for i = 1:size(exNumeric,1)
    for j = 1:size(exNumeric,2)
    fprintf(fOut, '%15.6e', exNumeric(i, j));
    end
    fprintf(fOut, '\n');
end
fclose(fOut);
