%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
%1) for normal-RHC and normal, long domain
xCoord= 5.014970e-02; % depend on centerline grid location
yMin = 0;
yMax = 0.5;
%%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%2) for normal-RHC and normal, cubic domain
yMin =  0;
yMax =  0.1;

grids = 5; % 5 10 15 20 25 30
% xCoord copied from line script;
% zCoord copied from line data.
if grids == 5
    xCoord = 5.059524e-02;
elseif grids == 10
    xCoord = 5.000000e-02;
elseif grids == 15
    xCoord = 5.020000e-02;
elseif grids == 20
    xCoord = 5.014970e-02;
elseif grids == 25
    xCoord = 5.000000e-02;
elseif grids == 30
    xCoord = 5.010000e-02;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yCoord= xCoord;
fileInfo = dir('couple_fluidplot_*.dat');
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    cfd_fluidplot_line(fileInfo(i).name, xCoord, yCoord, yMin, yMax);
end
