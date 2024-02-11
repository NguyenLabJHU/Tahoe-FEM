function cfd_fluidplot_line_sphericalwave(fileToRead1, xCoord, yCoord, yMin, yMax)

%newData1 = importdata(fileToRead1); % it is very slow
%numeric = newData1.data;
fid = fopen(fileToRead1);
text = fscanf(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 28);
numeric = fscanf(fid, '%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e', [21 inf]);
fclose(fid);
numeric = numeric.';
centerline = numeric(( numeric(:, 1) == xCoord & numeric(:, 2) == yCoord ), :);

colheaders = {'x' 'y' 'z' 'mach' 'density' 'momentumX' 'momentumY' 'momentumZ' ...
              'energy' 'velocityX' 'velocityY' 'velocityZ' 'pressure' 'temperature' ...
              'mask' 'penalFx' 'penalFy' 'penalFz' 'pressureFx' 'pressureFy' 'pressureFz'};
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
penalFz_c = evalin('caller', 'penalFz_c');
pressureFz_c = evalin('caller', 'pressureFz_c');
internalEnergy_c = energy_c./density_c-0.5*(velocityX_c.^2+velocityY_c.^2+velocityZ_c.^2);

path = pwd;
cell = strsplit(path, filesep);

screenWidth = 0.5; % work well for remote and local for dual screen
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1.0], 'visible', 'off');
set(0, 'CurrentFigure', fh);

scale = 0; % 0: xlimit set by Matlab; 1: xlimit specified
subplot(2,3,1);
%varX1 = pressure;
%varY1 = z;
varX2 = pressure_c;
varY2 = z_c;
%h=plot(varX1, varY1, 'r-', varX2, varY2, 'bo', 'LineWidth', 2);
h=plot(varX2, varY2, 'b-', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('pressure');
ylabel('radius');
ylim([yMin, yMax]);
if scale == 1
  xlim([0, 2.5e10]);
end

subplot(2,3,2);
%varX1 = density;
%varY1 = z;
varX2 = density_c;
varY2 = z_c;
%h=plot(varX1, varY1, 'r-', varX2, varY2, 'bo', 'LineWidth', 2);
h=plot(varX2, varY2, 'b-', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('density');
%ylabel('radius');
ylim([yMin, yMax]);
if scale == 1
  xlim([0, 2300]); %2000
end

subplot(2,3,3);
%varX1 = velocityZ;
%varY1 = z;
varX2 = velocityZ_c;
varY2 = z_c;
%h=plot(varX1, varY1, 'r-', varX2, varY2, 'bo', 'LineWidth', 2);
h=plot(varX2, varY2, 'b-', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('velocityZ');
%ylabel('radius');
ylim([yMin, yMax]);
if scale == 1
 xlim([-1e4, 1e4]); %6e3
end

subplot(2,3,4);
%varX1 = mach;
%varY1 = z;
varX2 = mach_c;
varY2 = z_c;
%h=plot(varX1, varY1, 'r-', varX2, varY2, 'b-', 'LineWidth', 2);
h=plot(varX2, varY2, 'b-', 'LineWidth', 2); % 'bo'
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('Mach');
ylabel('radius');
ylim([yMin, yMax]);
if scale == 1
  xlim([0, 8]); %5
end

subplot(2,3,5);
%varX1 = momentumZ;
%varY1 = z;
varX2 = momentumZ_c;
varY2 = z_c;
%h=plot(varX1, varY1, 'r-', varX2, varY2, 'bo', 'LineWidth', 2);
h=plot(varX2, varY2, 'b-', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('momentumZ');
%ylabel('radius');
ylim([yMin, yMax]);
if scale == 1
  xlim([-7e6, 7e6]); %3e6
end

subplot(2,3,6);
%varX1 = energy;
%varY1 = z;
varX2 = internalEnergy_c;
varY2 = z_c;
%h=plot(varX1, varY1, 'r-', varX2, varY2, 'bo', 'LineWidth', 2);
h=plot(varX2, varY2, 'b-', 'LineWidth', 2);
%set(h(1),'LineWidth',1);
%set(h(2),'LineWidth',2);
title('internal energy');
%ylabel('radius');
ylim([yMin, yMax]);
if scale == 1
  xlim([0, 3.5e7]);
end

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.line.png'), 'png');
%print(fh, '-r600', '-dpng', strcat(fileToRead1, '.png'));
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '.line.png'), options);

%curve
%dlmwrite('point.dat', curve, 'delimiter', '\t', 'precision', 6);
fOut=fopen(strcat(fileToRead1, '.line'), 'w');
fprintf(fOut, '%15s%15s%15s%15s%15s%15s%15s\n', 'location', 'pressure', 'density', 'velocityZ', 'Mach', 'momentumZ', 'internalEnergy');
for i = 1:size(centerline, 1)
    fprintf(fOut, '%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e\n', z_c(i), pressure_c(i), density_c(i), velocityZ_c(i), mach_c(i), momentumZ_c(i), internalEnergy_c(i));
end
fclose(fOut);

close(fh);
clear;
