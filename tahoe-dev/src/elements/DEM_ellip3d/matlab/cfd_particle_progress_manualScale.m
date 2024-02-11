function cfd_particle_progress(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % work well for remote and local for dual screen

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

accruedTime = evalin('caller', 'accruedTime');
penalFz     = evalin('caller', 'penalFz');
pressureFz  = evalin('caller', 'pressureFz');
internalFz  = evalin('caller', 'internalFz');
totalFz     = penalFz + pressureFz + internalFz;
accelZ      = evalin('caller', 'accelZ');
velocZ      = evalin('caller', 'velocZ');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

subplot(1,2,1)
forceH = plot(accruedTime, penalFz, 'r--o', accruedTime, pressureFz, 'm--', accruedTime, internalFz, 'b-x', accruedTime, totalFz, 'k-p','LineWidth', 2); 
title('drag force');
xlabel('time (s)');
ylabel('force (N)');
legend('viscousFz', 'pressureFz', 'internalFz', 'totalFz'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);
ylim([0, 800]);
%ylim([0, 45]);

subplot(1,2,2)
[AX,H1,H2] = plotyy(accruedTime, accelZ, accruedTime, velocZ);
set(H1, 'Color', 'b', 'LineStyle', '--', 'Marker', 'o', 'LineWidth', 2);
set(H2, 'Color', 'r', 'LineStyle', '-',  'Marker', 'x', 'LineWidth', 2);
set(AX(1),'YColor','b');
set(AX(2),'YColor','r');
title('kinematics');
xlabel('time (s)'); 
set(get(AX(1), 'Ylabel'), 'String', 'accelZ (m/s^2)');
set(get(AX(2), 'Ylabel'), 'String', 'velocZ (m/s)');
%set(AX(1), 'xlim', [0, 3.5e-4]);
%set(AX(2), 'xlim', [0, 0.4e-4]);

set(AX(1), 'ylim', [0, 1.0e+6]);
set(AX(1), 'ytick', [0:1.0e+5:1.0e+6]);
set(AX(2), 'ylim', [0, 25]);
set(AX(2), 'ytick', [0:5:25]);

%{
set(AX(1), 'ylim', [0, 7.5e+4]);
set(AX(1), 'ytick', [0:1.5e+4:7.5e+4]);
set(AX(2), 'ylim', [0, 12]);
set(AX(2), 'ytick', [0:1:12]);
legend('accelZ', 'velocZ'); % , 'location', 'best'
%}

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '.png'), options);
