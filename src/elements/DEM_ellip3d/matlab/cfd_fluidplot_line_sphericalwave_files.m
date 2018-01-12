charge = 2;
if charge == 1
%for spherical wave shock from explosive charge r=0.5 cm
 xCoord = 1.25e-3;
 yCoord = 1.25e-3; 
 yMin   = -2.1E-1;
 yMax   =  2.1E-1;
else if charge == 2
%for spherical wave shock from explosive charge r=2.5 cm
 xCoord = 2.000681e-01;
 yCoord = 2.000681e-01; 
 yMin   = 1.0E-1;
 yMax   = 6.5E-1;
end
end

fileInfo = dir('couple_fluidplot_*.dat');
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    cfd_fluidplot_line_sphericalwave(fileInfo(i).name, xCoord, yCoord, yMin, yMax);
end
