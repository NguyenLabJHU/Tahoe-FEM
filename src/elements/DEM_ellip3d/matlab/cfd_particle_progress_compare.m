%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % work well for remote and local for dual screen

data1 = importdata('particle_progress');
for i = 1:size(data1.colheaders, 2)
    assignin('base', genvarname(data1.colheaders{i}), data1.data(:,i));
end

data2 = importdata('particle_progress_1.00');
for i = 1:size(data2.colheaders, 2)
    assignin('base', strcat(genvarname(data2.colheaders{i}), '_c'), data2.data(:,i));
end

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

forceH = plot(accruedTime, penalFz, 'r-o', accruedTime_c, penalFz_c, 'g-x', ...
              accruedTime, pressureFz, 'b-s', accruedTime_c, pressureFz_c, 'm-*','LineWidth', 2);
title('Penalization drag force and pressure acting on a stationary particle');
xlabel('time (s)');
ylabel('forces (N)');
legend('viscousFz-mass penalization', 'viscousFz', 'pressureFz-mass penalization', 'pressureFz', 'location', 'East'); %, 'location', 'best'
grid on
box on
axis square
%axis([1.5e-5, 4.0e-5, ylim]);
%xlim([0, 3.5e-4]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, 'particle_progress-masspenl.png', 'png');
options.Format = 'png';
hgexport(fh, 'particle_progress-masspenl.png', options);
