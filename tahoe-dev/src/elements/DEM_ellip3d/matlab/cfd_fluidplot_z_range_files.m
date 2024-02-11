xGrid=19; % grid number in x direction
yGrid=19; % grid number in y direction
z1Coord=-7.5e-2;
z2Coord= 7.5e-2;

fileInfo = dir('couple_fluidplot_*.dat');
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    cfd_fluidplot_z_range(fileInfo(i).name, xGrid, yGrid, z1Coord, z2Coord);
end
