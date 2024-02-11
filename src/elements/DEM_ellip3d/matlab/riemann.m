function riemann(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % work well for remote and local for dual screen

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

position = evalin('caller', 'position');
density  = evalin('caller', 'density');
velocity = evalin('caller', 'velocity');
momentum = evalin('caller', 'momentum');
pressure = evalin('caller', 'pressure');
int_energy = evalin('caller', 'int_energy');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'off');

subplot(1,5,1)
plot(density, position, 'r-', 'LineWidth', 2);
title('density');
ylabel('position');
ylim([0, 1]);

subplot(1,5,2)
plot(velocity, position, 'r-', 'LineWidth', 2);
title('velocity');
ylim([0, 1]);

subplot(1,5,3)
plot(momentum, position, 'r-', 'LineWidth', 2);
title('momentum');
ylim([0, 1]);

subplot(1,5,4)
plot(pressure, position, 'r-', 'LineWidth', 2);
title('pressure');
ylim([0, 1]);

subplot(1,5,5)
plot(int_energy, position, 'r-', 'LineWidth', 2);
title('internal energy');
ylim([0, 1]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
saveas(fh, strcat(fileToRead1, '.png'), 'png');

close(fh);
clear;
