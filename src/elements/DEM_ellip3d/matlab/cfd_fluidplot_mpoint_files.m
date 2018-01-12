screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%1) for shock tube experiments of UT Dallas
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
%{
%2) for spherical shock wave w/o particles, used for calibration.
%   not comparable to cases with particles due to different spatial positioning.
timeTotal = 3.0e-04;
timeSteps = 100;
timeSnap = timeTotal / timeSteps;
xCoord = 1.25e-3;
yCoord = 1.25e-3; 
zCoord = [1.125000e-02; 2.125000e-02; 3.875000e-02; 4.875000e-02; 6.875000e-02; 1.012500e-01; 1.512500e-01; 2.012500e-01];
suffix = {'_1_r1.1e-2'; '_2_r2.1e-2'; '_3_r3.9e-2'; '_4_r4.9e-2'; '_5_r6.9e-2'; '_6_r1.0e-1'; '_7_r1.5e-1'; '_8_r2.0e-1'}; % cell array
%}
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
timeTotal = 3.0e-04;
timeSteps = 100;
timeSnap = timeTotal / timeSteps;
xCoord = 1.007463e-01;
yCoord = 1.007463e-01;
% coordinates x y z are estimated in terms of tracking particles such as V1, V2, V3 ... or
% some heights above soil surface, however they must be read from couple_fluidplot_000.dat 
% by Matlab script extract_fluidpoints.m for accurate numbers.
dob = 5.1; %5.1, 7.6
if dob == 5.1
  zCoord = [1.395522e-01; 1.455224e-01; 1.708955e-01; 1.962687e-01; 2.216418e-01; 2.470149e-01];
  suffix = {'_0_z1.39e-1_charge_center'; '_1_z1.46e-1_r0.60cm_charge_edge'; '_2_z1.71e-1_r3.13cm'; ... 
            '_3_z1.96e-1_r5.67cm_ground';'_4_z2.21e-01_r8.21cm_above_2.54cm';'_5_z2.47e-1_r10.75cm_above_5.07cm'}; % cell array
elseif dob == 7.6
  zCoord = [1.141791e-01; 1.201493e-01; 1.455224e-01; 1.708955e-01; 1.962687e-01];
  suffix = {'_0_z1.14e-1_charge_center'; '_1_z1.20e-1_r0.61cm_charge_edge'; '_2_z1.45e-1_r3.2cm'; '_3_z1.71e-1_r5.7cm'; '_4_z1.96e-1_r8.2cm_ground';}; % cell array
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
%6) for spherical wave impacting thousands of particle from inside a cavity
% 18000 particles
timeTotal = 10e-3; %0.3e-3;
timeSteps = 100;
timeSnap  = timeTotal / timeSteps;
xCoord = 1.505618e-01;
yCoord = 1.505618e-01;
% coordinates x y z are estimated in terms of tracking particles such as V1, V2, V3 ... or
% some heights above soil surface, however they must be read from couple_fluidplot_000.dat 
% by Matlab script extract_fluidpoints.m for accurate numbers.
dob = 10.1; %5.1, 7.6 10.1 cm
if dob == 5.1
  zCoord = [1.801498e-01; 2.041199e-01; 2.310861e-01; 2.580524e-01; ... 
            2.835206e-01; 3.044944e-01; 3.254682e-01; ...
            2.745318e-01; ...
            3.359551e-01; 3.853933e-01; 4.348315e-01; 4.857678e-01];
  suffix = {'_V7_z1.80e-1_r-9.4cm'; '_V6_z2.04e-1_r-7.0cm'; '_V5_z2.31e-1_r-4.3cm'; '_V4_z2.58e-1_r-1.65cm'; ... 
            '_V3_z2.84e-01_r0.9cm'; '_V2_z3.04e-1_r3.0cm';  '_V1_z3.25e-01_r5.1cm'; ...
            '_D1_z2.745318e-01_charge_center'; ...
            '_A1_z3.36e-1_r6.1cm_ground'; '_A2_z3.85e-01_r11.1cm_above_5cm'; '_A3_z4.35e-1_r16.0cm_above_10cm'; '_A4_z4.86e-1_r21.1cm_above_15cm'}; % cell array
elseif dob == 7.6
  zCoord = [1.801498e-01; 2.041199e-01; 2.310861e-01; 2.580524e-01; ... 
            2.835206e-01; 3.044944e-01; 3.254682e-01; ...
            2.490637e-01; ...
            3.359551e-01; 3.853933e-01; 4.348315e-01; 4.857678e-01];
  suffix = {'_V7_z1.80e-1_r-6.9cm'; '_V6_z2.04e-1_r-4.5cm'; '_V5_z2.31e-1_r-1.8cm'; '_V4_z2.58e-1_r0.9cm'; ... 
            '_V3_z2.84e-01_r3.5cm'; '_V2_z3.04e-1_r5.5cm';  '_V1_z3.25e-01_r7.6cm'; ...
            '_D2_z2.490637e-01_charge_center'; ...
            '_A1_z3.36e-1_r8.7cm_ground'; '_A2_z3.85e-01_r13.7cm_above_5cm'; '_A3_z4.35e-1_r18.6cm_above_10cm'; '_A4_z4.86e-1_r23.7cm_above_15cm'}; % cell array
elseif dob == 10.1
  zCoord = [1.801498e-01; 2.041199e-01; 2.310861e-01; 2.580524e-01; ... 
            2.835206e-01; 3.044944e-01; 3.254682e-01; ...
            2.235955e-01; ...
            3.359551e-01; 3.853933e-01; 4.348315e-01; 4.857678e-01];
  suffix = {'_V7_z1.80e-1_r-4.3cm'; '_V6_z2.04e-1_r-2.0cm'; '_V5_z2.31e-1_r0.8cm'; '_V4_z2.58e-1_r3.5cm'; ... 
            '_V3_z2.84e-01_r6.0cm'; '_V2_z3.04e-1_r8.1cm';  '_V1_z3.25e-01_r10.2cm'; ...
            '_D3_2.235955e-01_charge_center'; ...
            '_A1_z3.36e-1_r11.2cm_ground'; '_A2_z3.85e-01_r16.2cm_above_5cm'; '_A3_z4.35e-1_r21.1cm_above_10cm'; '_A4_z4.86e-1_r26.2cm_above_15cm'}; % cell array
end
%%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%7) for spherical wave impacting thousands of particle from inside a cavity
% 48000 particles
timeTotal = 1.5e-3; %3.0e-4;
timeSteps = 100;
timeSnap  = timeTotal / timeSteps;
% coordinates x y z are estimated in terms of tracking particles such as V1, V2, V3 ... or
% some heights above soil surface, however they must be read from couple_fluidplot_000.dat 
% by Matlab script extract_fluidpoints.m for accurate numbers.
xCoord = 2.000681e-01;
yCoord = 2.000681e-01;
dob = 29; %10, 20, 29
if dob == 10
  zCoord = [1.711853e-01; 2.626022e-01; 3.465259e-01; 3.839918e-01; 4.319482e-01; ... 
            3.150545e-01; ...
            4.499319e-01; 5.008856e-01; 5.503406e-01; 6.012943e-01; 6.492507e-01];
  suffix = {'_V5_z1.711853e-01_r-14.4cm'; '_V4_z2.626022e-01_r-5.3cm'; '_V3_z3.465259e-01_r3.2cm'; '_V2_z3.839918e-01_r6.9cm'; '_V1_z4.319482e-01_r11.7cm'; ...
            '_D1_z3.150545e-01_charge_center'; ...
            '_A1_z4.499319e-01_ground'; '_A2_z5.008856e-01_above_5cm'; '_A3_z5.503406e-01_above_10cm'; '_A4_z6.012943e-01_above_15cm'; '_A5_z6.492507e-01_above_20cm'}; % cell array
elseif dob == 20
  zCoord = [1.711853e-01; 2.626022e-01; 3.465259e-01; 3.839918e-01; 4.319482e-01; ... 
            2.146458e-01; ...
            4.499319e-01; 5.008856e-01; 5.503406e-01; 6.012943e-01; 6.492507e-01];
  suffix = {'_V5_z1.711853e-01_r-4.4cm'; '_V4_z2.626022e-01_r4.8cm'; '_V3_z3.465259e-01_r13.2cm'; '_V2_z3.839918e-01_r16.9cm'; '_V1_z4.319482e-01_r21.7cm'; ...
            '_D2_z2.146458e-01_charge_center'; ...
            '_A1_z4.499319e-01_ground'; '_A2_z5.008856e-01_above_5cm'; '_A3_z5.503406e-01_above_10cm'; '_A4_z6.012943e-01_above_15cm'; '_A5_z6.492507e-01_above_20cm'}; % cell array
elseif dob == 29
  zCoord = [1.711853e-01; 2.626022e-01; 3.465259e-01; 3.839918e-01; 4.319482e-01; ... 
            1.247275e-01; ...
            4.499319e-01; 5.008856e-01; 5.503406e-01; 6.012943e-01; 6.492507e-01];
  suffix = {'_V5_z1.711853e-01_r4.7cm'; '_V4_z2.626022e-01_r13.8cm'; '_V3_z3.465259e-01_r22.2cm'; '_V2_z3.839918e-01_r25.9cm'; '_V1_z4.319482e-01_r30.7cm'; ...
            '_D3_z1.247275e-01_charge_center'; ...
            '_A1_z4.499319e-01_ground'; '_A2_z5.008856e-01_above_5cm'; '_A3_z5.503406e-01_above_10cm'; '_A4_z6.012943e-01_above_15cm'; '_A5_z6.492507e-01_above_20cm'}; % cell array
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileInfo = dir('couple_fluidplot_*.dat');
curve = zeros(size(zCoord,1), size(fileInfo, 1), 7);
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    out=zeros(size(zCoord,1), 6); 
    [out] = cfd_fluidplot_mpoint(fileInfo(i).name, xCoord, yCoord, zCoord);
    for j = 1:size(zCoord,1)
        curve(j, i, 1) = i * timeSnap;
        curve(j, i, 2:7) = out(j, 1:6);
    end
end

%curve
%dlmwrite('point.dat', curve, 'delimiter', '\t', 'precision', 6);

for k = 1:size(zCoord,1)
fname= strcat(strcat('point', char(suffix(k))), '.dat');

fOut=fopen(fname, 'w');
fprintf(fOut, '%15s%15s%15s%15s%15s%15s%15s\n', 'time', 'pressure', 'density', 'velocityZ', 'Mach_number', 'total_energy', 'temperature');
for i = 1:size(curve,2)
    for j = 1:size(curve,3)
    fprintf(fOut, '%15.6e', curve(k, i, j));
    end
    fprintf(fOut, '\n');
end
fclose(fOut);

%moved to top
%screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1.0], 'visible', 'off');
set(0, 'CurrentFigure', fh); % 0 is monitor handle

plottype = 0; % 0: plot 6 variables; 1: plot 3 variables
if plottype == 0

subplot(2,3,1);
h=plot(curve(k, :, 1), curve(k, :, 2), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('pressure');
xlim([0, timeTotal]);

subplot(2,3,2);
h = plot(curve(k, :, 1), curve(k, :, 3), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('density');
xlim([0, timeTotal]);

subplot(2,3,3);
h = plot(curve(k, :, 1), curve(k, :, 4), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('velocityZ');
xlim([0, timeTotal]);

subplot(2,3,4);
h = plot(curve(k, :, 1), curve(k, :, 5), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('Mach Number');
xlim([0, timeTotal]);

subplot(2,3,5);
h = plot(curve(k, :, 1), curve(k, :, 6), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('total energy');
xlim([0, timeTotal]);

subplot(2,3,6);
h = plot(curve(k, :, 1), curve(k, :, 7), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('temperature');
xlim([0, timeTotal]);

elseif plottype == 1
subplot(1,3,1);
h = plot(curve(k, :, 1), curve(k, :, 2), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('pressure');
xlim([0, timeTotal]);

subplot(1,3,2);
h = plot(curve(k, :, 1), curve(k, :, 3), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('density');
xlim([0, timeTotal]);

subplot(1,3,3);
h = plot(curve(k, :, 1), curve(k, :, 4), 'b-o', 'LineWidth', 2);
xlabel('time');
ylabel('velocityZ');
xlim([0, timeTotal]);

suffix(k) = strcat(char(suffix(k)), '_dup');
end

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(strcat('point', suffix(k)), '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(strcat('point', char(suffix(k))), '.png'), options);

end
