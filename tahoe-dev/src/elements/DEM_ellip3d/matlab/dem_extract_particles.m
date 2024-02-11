function dem_extract_particles(fileToRead1)

x1Coord= 2.45e-1;
x2Coord= 2.55e-1;
%y1Coord= 9.5e-3;
%y2Coord= 10e-3;
%z1Coord= 12e-3;
%z2Coord= 13e-3;

fIn = fopen(fileToRead1);
text  = fscanf(fIn, '%s', 30);
numeric = fscanf(fIn, '%e', [29 inf]);
fclose(fIn);
numeric = numeric.';
exNumeric = numeric(( numeric(:, 6) >= x1Coord & numeric(:, 6) <= x2Coord ), :);

fOut=fopen(strcat('ex_', fileToRead1), 'w');

fprintf(fOut, '%15d\n', size(exNumeric,1));
fprintf(fOut, '%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n', ...
              'id', 'type', 'radius_a', 'radius_b', 'radius_c', 'position_x', 'position_y', 'position_z', ...
              'axis_a_x', 'axis_a_y', 'axis_a_z', 'axis_b_x', 'axis_b_y', 'axis_b_z', 'axis_c_x', 'axis_c_y', 'axis_c_z', ...
              'velocity_x', 'velocity_y', 'velocity_z', 'omga_x', 'omga_y', 'omga_z', ...
              'force_x', 'force_y', 'force_z', 'moment_x', 'moment_y', 'moment_z');

for i = 1:size(exNumeric,1)
    fprintf(fOut, '%15.6d%15.6d', exNumeric(i, 1), exNumeric(i, 2));
    for j = 3:size(exNumeric,2)
    fprintf(fOut, '%15.6e', exNumeric(i, j));
    end
    fprintf(fOut, '\n');
end
fclose(fOut);
