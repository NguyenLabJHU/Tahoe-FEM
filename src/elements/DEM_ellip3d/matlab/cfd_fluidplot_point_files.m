%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%1) for shock tube experiments
timeSnap = 20000 * 5.0E-7 / 100; %7500 * 1.0E-8 / 100;
xCoord =  0.000000e+00;
yCoord = -1.042105e-02;
zCoord = -6.606551e-02; % sensor Ch1
%zCoord = -2.968671e-03; % sensor Ch2, 30:-2.970000e-03, 40:-2.969397e-03, 50:-2.968978e-03, "60":-2.968671e-03, 80:-2.968250e-03
%zCoord = 5.364873e-03; % sensor Ch3
%zCoord = 1.250791e-02; % sensor Ch4p
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
%2) for spherical shock wave w/o particles
timeTotal = 2.0e-04;
timeSteps = 100;
zCoordChoice = 3; % 1 2 3 4 5 6 from smaller to larger radius

timeSnap = timeTotal / timeSteps;
xCoord = 1.25e-3;
yCoord = 1.25e-3; 

if zCoordChoice == 1
    zCoord = 1.125000e-02;
    suffix = '_1_r1.1e-2';
elseif zCoordChoice == 2
    zCoord = 2.125000e-02;
    suffix = '_2_r2.1e-2';
elseif zCoordChoice == 3
    zCoord = 4.875000e-02;
    suffix = '_3_r4.9e-2';
elseif zCoordChoice == 4
    zCoord = 1.012500e-01;
    suffix = '_4_r1.0e-1';
elseif zCoordChoice == 5
    zCoord = 1.512500e-01;
    suffix = '_5_r1.5e-1';
elseif zCoordChoice == 6
    zCoord = 2.012500e-01;
    suffix = '_6_r2.0e-1';
end
%%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%3) for normal-RHC and normal for long domain
% changeable parameters:
timeSnap = 1600 * 5.0E-7 / 100;
zCoordChoice = 3; % 1, 2, 3 from lower to higher location in z direction.

% fixed parameters:
xCoord = 5.014970e-02;
yCoord = 5.014970e-02;
if zCoordChoice == 1
    zCoord = 5.239521e-03;
    suffix = '_1_z5.2e-3';
elseif zCoordChoice == 2
    zCoord = 1.010479e-01;
    suffix = '_2_z1.0e-1';
elseif zCoordChoice == 3
    zCoord = 4.004491e-01;
    suffix = '_3_z4.0e-1';
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%4) for normal-RHC and normal for cubic domain
% changeable parameters:
timeTotal = 2.5e-5;
timeSteps = 25;
grids = 10; % 5 10 15 20 25 30

timeSnap  = timeTotal / timeSteps;
% xCoord copied from line script,
% zCoord copied from line data.
if grids == 5
    xCoord = 5.059524e-02;
    zCoord = 5.297619e-02;
elseif grids == 10
    xCoord = 5.000000e-02;
    zCoord = 5.239521e-02;
elseif grids == 15
    xCoord = 5.020000e-02;
    zCoord = 5.220000e-02;
elseif grids == 20
    xCoord = 5.014970e-02;
    zCoord = 5.224551e-02;
elseif grids == 25
    xCoord = 5.000000e-02;
    zCoord = 5.191847e-02;
elseif grids == 30
    xCoord = 5.010000e-02;
    zCoord = 5.210000e-02;
end

yCoord = xCoord;
suffix = '_stagnation';
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%5) for spherical wave impacting thousands of particle from inside a cavity 
% 3000 particles
dob = 7.9; %5.1, 7.9
timeTotal = 4.0e-4;
timeSteps = 100;

timeSnap  = timeTotal / timeSteps;
xCoord = 1.007463e-01;
yCoord = 1.007463e-01;
if dob == 5.1
    zCoord = 1.500000e-01; % edge of the explosive charge
elseif dob == 7.9
    zCoord = 1.216418e-01; % edge of the explosive charge
end
suffix = '_charge_edge';
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%6) for spherical wave impacting thousands of particle from inside a cavity
% 14000 particles
timeTotal = 4.0e-4;
timeSteps = 100;
timeSnap  = timeTotal / timeSteps;
xCoord = 1.011236e-01;
yCoord = 1.011236e-01;
dob = 15; %7.9, 15
if dob == 7.9
    zCoord = 0; % center edge of the explosive charge
elseif dob == 15
    zCoord = 1.216418e-01; % center edge of the explosive charge
end
suffix = '_charge_edge';
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileInfo = dir('couple_fluidplot_*.dat');
curve=zeros(size(fileInfo, 1), 7);
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    curve(i, 1) = i * timeSnap;
    [curve(i, 2), curve(i, 3), curve(i, 4), curve(i, 5), curve(i, 6), curve(i, 7)] = cfd_fluidplot_point(fileInfo(i).name, xCoord, yCoord, zCoord);
end
%curve
%dlmwrite('point.dat', curve, 'delimiter', '\t', 'precision', 6);
fOut=fopen(strcat(strcat('point', suffix), '.dat'), 'w');
fprintf(fOut, '%15s%15s%15s%15s%15s%15s%15s\n', 'time', 'pressure', 'density', 'velocityZ', 'Mach_number', 'total_energy', 'temperature');
for i = 1:size(curve,1)
    for j = 1:size(curve,2)
    fprintf(fOut, '%15.6e', curve(i, j));
    end
    fprintf(fOut, '\n');
end
fclose(fOut);

screenWidth = 0.5; % work well for remote and local for dual screen
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1.0], 'visible', 'off');
set(0, 'CurrentFigure', fh); % 0 is monitor handle

plottype = 0; % 0: plot 6 variables; 1: plot 3 variables
if plottype == 0

subplot(2,3,1);
h=plot(curve(:, 1), curve(:, 2), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('pressure');
xlim([0, timeTotal]);

subplot(2,3,2);
h = plot(curve(:, 1), curve(:, 3), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('density');
xlim([0, timeTotal]);

subplot(2,3,3);
h = plot(curve(:, 1), curve(:, 4), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('velocityZ');
xlim([0, timeTotal]);

subplot(2,3,4);
h = plot(curve(:, 1), curve(:, 5), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('Mach Number');
xlim([0, timeTotal]);

subplot(2,3,5);
h = plot(curve(:, 1), curve(:, 6), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('total energy');
xlim([0, timeTotal]);

subplot(2,3,6);
h = plot(curve(:, 1), curve(:, 7), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('temperature');
xlim([0, timeTotal]);

elseif plottype == 1
subplot(1,3,1);
h = plot(curve(:, 1), curve(:, 2), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('pressure');
xlim([0, timeTotal]);

subplot(1,3,2);
h = plot(curve(:, 1), curve(:, 3), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('density');
xlim([0, timeTotal]);

subplot(1,3,3);
h = plot(curve(:, 1), curve(:, 4), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('velocityZ');
xlim([0, timeTotal]);

suffix = strcat(suffix, '_dup');
end

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(strcat('point', suffix), '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(strcat('point', suffix), '.png'), options);
