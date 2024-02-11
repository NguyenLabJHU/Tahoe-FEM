function cfd_particle_progress_allforce(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % work well for remote and local for dual screen

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

accruedTime = evalin('caller', 'accruedTime');
penalFx     = evalin('caller', 'penalFx');
penalFy     = evalin('caller', 'penalFy');
penalFz     = evalin('caller', 'penalFz');
pressureFx  = evalin('caller', 'pressureFx');
pressureFy  = evalin('caller', 'pressureFy');
pressureFz  = evalin('caller', 'pressureFz');
penalMx     = evalin('caller', 'penalMx');
penalMy     = evalin('caller', 'penalMy');
penalMz     = evalin('caller', 'penalMz');
pressureMx  = evalin('caller', 'pressureMx');
pressureMy  = evalin('caller', 'pressureMy');
pressureMz  = evalin('caller', 'pressureMz');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

subplot(1,2,1)
forceH = plot(accruedTime, penalFx, 'r-', accruedTime, penalFy, 'b-', accruedTime, penalFz, 'g-', ...
              accruedTime, pressureFx, 'm--', accruedTime, pressureFy, 'k--', accruedTime, pressureFz, 'c--', 'LineWidth', 2); % o x
title('force');
xlabel('time (s)');
ylabel('force (N)');
legend('viscousFx', 'viscousFy', 'viscousFz', 'pressureFx', 'pressureFy', 'pressureFz'); %, 'location', 'best'
%axis([0.2e-4, 0.3e-4, ylim]);

subplot(1,2,2)
forceH = plot(accruedTime, penalMx, 'r-', accruedTime, penalMy, 'b-', accruedTime, penalMz, 'g-', ...
              accruedTime, pressureMx, 'm--', accruedTime, pressureMy, 'k--', accruedTime, pressureMz, 'c--', 'LineWidth', 2); % o x
title('moment');
xlabel('time (s)');
ylabel('moment (N.m)');
legend('viscousMx', 'viscousMy', 'viscousMz', 'pressureMx', 'pressureMy', 'pressureMz'); %, 'location', 'best'
%axis([0.2e-4, 0.3e-4, ylim]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '_allforce.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_allforce.png'), options);
