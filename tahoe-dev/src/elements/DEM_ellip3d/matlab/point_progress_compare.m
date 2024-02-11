%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % work well for remote and local for dual screen

data1 = importdata('point_progress_ch1.dat');
for i = 1:size(data1.colheaders, 2)
    assignin('base', genvarname(data1.colheaders{i}), data1.data(:,i));
end

data2 = importdata('point_progress_ch2.dat');
for i = 1:size(data2.colheaders, 2)
    assignin('base', strcat(genvarname(data2.colheaders{i}), '_c'), data2.data(:,i));
end

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

subplot(1,3,1)
plot(time, pressure, 'r-o', time_c, pressure_c, 'b-x', 'LineWidth', 2);
title('');
xlabel('time');
ylabel('pressure');
legend('pressure-Ch1', 'pressure-Ch2', 'location', 'best'); %, 'location', 'best'
grid on
box on
%axis([0.1e-4, 8.0e-5, ylim]);
%xlim([0, 3.5e-4]);

subplot(1,3,2)
plot(time, density, 'r-o', time_c, density_c, 'b-x', 'LineWidth', 2);
title('');
xlabel('time');
ylabel('density');
legend('density-Ch1', 'density-Ch2', 'location', 'best'); %, 'location', 'best'
grid on
box on
%axis([0.1e-4, 8.0e-5, ylim]);
%xlim([0, 3.5e-4]);

subplot(1,3,3)
plot(time, velocityZ, 'r-o', time_c, velocityZ_c, 'b-x', 'LineWidth', 2);
title('');
xlabel('time');
ylabel('velocityZ');
legend('velocityZ-Ch1', 'velocityZ-Ch2', 'location', 'best'); %, 'location', 'best'
grid on
box on
%axis([0.1e-4, 8.0e-5, ylim]);
%xlim([0, 3.5e-4]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
saveas(fh, 'point_progress-shocktube-noptcl.png', 'png');
