function cfd_fluidplot_z_range(fileToRead1, xGrid, yGrid, z1Coord, z2Coord)

fIn = fopen(fileToRead1);
text  = fscanf(fIn, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 28);
numeric = fscanf(fIn, '%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e', [21 inf]);
fclose(fIn);
numeric = numeric.';
exNumeric = numeric(( numeric(:, 3) >= z1Coord & numeric(:, 3) <= z2Coord ), :);

fOut=fopen(strcat('ex_', fileToRead1), 'w');
fprintf(fOut, '%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s%15s\n', ...
              'VARIABLES = "x"', '"y"', '"z"', '"mach"', '"density"', '"momentumX"', '"momentumY"', '"momentumZ"', ...
              '"energy"', '"velocityX"', '"velocityY"', '"velocityZ"', '"pressure"', '"temperature"', ...
              '"mask"', '"penalFx"', '"penalFy"', '"penalFz"', '"pressureFx"', '"pressureFy"', '"pressureFz"');
fprintf(fOut, 'ZONE I=%d, J=%d, K=%d, DATAPACKING=POINT\n', xGrid, yGrid, size(exNumeric, 1)/(xGrid*yGrid));
for i = 1:size(exNumeric,1)
    for j = 1:size(exNumeric,2)
    fprintf(fOut, '%15.6e', exNumeric(i, j));
    end
    fprintf(fOut, '\n');
end
fclose(fOut);


