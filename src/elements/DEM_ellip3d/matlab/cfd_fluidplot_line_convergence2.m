function cfd_fluidplot_line_convergence()
% for cubic domain
grids = 5; % 5 10 15 20 25 30
zoomin = 0; % 0

fileToRead1='couple_fluidplot_020.dat';
if zoomin == 0
  zMin = 0;
  zMax = 0.1;
elseif zoomin == 1
  zMin = 0.045;
  zMax = 0.065;
end

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
yCoord = xCoord;
% end for cubic domain

%newData1 = importdata(fileToRead1); % it is very slow
%numeric = newData1.data;
fid = fopen(fileToRead1);
text = fscanf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 31);
numeric = fscanf(fid, '%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e', [24 inf]);
fclose(fid);
numeric = numeric.';
centerline = numeric((numeric(:, 1) == xCoord & numeric(:, 2) == yCoord & numeric(:, 3) >= zMin & numeric(:, 3)<= zMax), :);

colheaders = {'x' 'y' 'z' 'mach' 'density' 'momentumX' 'momentumY' 'momentumZ' ...
              'energy' 'velocityX' 'velocityY' 'velocityZ' 'pressure' 'temperature' ...
              'mask' 'viscPenalFx' 'viscPenalFy' 'viscPenalFz' 'presPenalFx' 'presPenalFy' 'presPenalFz' ...
              'pressureFx' 'pressureFy' 'pressureFz'};
for i = 1:size(numeric, 2)
    %assignin('caller', genvarname(colheaders{i}), numeric(:,i));
    assignin('caller', strcat(genvarname(colheaders{i}), '_c'), centerline(:,i));
end

%z = evalin('caller', 'z');
%energy = evalin('caller', 'energy');
%mach = evalin('caller', 'mach');
%density = evalin('caller', 'density');
%velocityZ  = evalin('caller', 'velocityZ');
%pressure = evalin('caller', 'pressure');
%temperature = evalin('caller', 'temperature');
%penalFz = evalin('caller', 'penalFz');
%pressureFz = evalin('caller', 'pressureFz');

z_c = evalin('caller', 'z_c');
energy_c = evalin('caller', 'energy_c');
mach_c = evalin('caller', 'mach_c');
density_c = evalin('caller', 'density_c');
velocityX_c  = evalin('caller', 'velocityX_c');
velocityY_c  = evalin('caller', 'velocityY_c');
velocityZ_c  = evalin('caller', 'velocityZ_c');
momentumZ_c  = evalin('caller', 'momentumZ_c');
pressure_c = evalin('caller', 'pressure_c');
temperature_c = evalin('caller', 'temperature_c');
penalFz_c = evalin('caller', 'viscPenalFz_c');
pressureFz_c = evalin('caller', 'pressureFz_c');
mask_c = evalin('caller', 'mask_c');

path = pwd;
cell = strsplit(path, filesep);

screenWidth = 0.5; % work well for remote and local for dual screen
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1.0]);%, 'visible', 'off');
%set(0, 'CurrentFigure', fh);

%{
subplot(1,2,1);
varX1 = mask_c * max(density_c);
%varY1 = z;
varX2 = density_c;
varY2 = z_c;
h=plot(varX2, varY2, 'b-o', varX1, varY2, 'r-*', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('density');
ylabel('position');
ylim([zMin, zMax]);
legend('density', 'mask');
%}


subplot(1,2,1);
varX1 = mask_c * max(pressure_c);
%varY1 = z;
varX2 = pressure_c;
varY2 = z_c;
h=plot(varX2, varY2, 'b-o', varX1, varY2, 'r-*', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('pressure');
ylabel('position');
ylim([zMin, zMax]);
legend('pressure', 'mask');
%xlim([0, 1.0e+7]);

subplot(1,2,2);
varX1 = mask_c * max(mach_c);
%varY1 = z;
varX2 = mach_c;
varY2 = z_c;
h=plot(varX2, varY2, 'b-o', varX1, varY2, 'r-*', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('Mach number');
ylabel('position');
ylim([zMin, zMax]);
legend('Mach number', 'mask');

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.line.png'), 'png');
%print(fh, '-r600', '-dpng', strcat(fileToRead1, '.png'));
options.Format = 'png';
if zoomin == 0
  hgexport(fh, strcat(fileToRead1, '.line.png'), options);
elseif zoomin == 1
  hgexport(fh, strcat(fileToRead1, '.line.zoomin.png'), options);
end

%curve
%dlmwrite('point.dat', curve, 'delimiter', '\t', 'precision', 6);
%%{
if zoomin == 0
  fOut=fopen(strcat(fileToRead1, '.line'), 'w');
elseif zoomin == 1
  fOut=fopen(strcat(fileToRead1, '.line.zoomin'), 'w');
end
fprintf(fOut, '%15s%15s%15s\n', 'location', 'Mach_number', 'pressure');
for i = 1:size(centerline, 1)
    fprintf(fOut, '%15.6e%15.6e%15.6e\n', z_c(i), mach_c(i), pressure_c(i));
end
fclose(fOut);
%%}
%close(fh);
clear;
