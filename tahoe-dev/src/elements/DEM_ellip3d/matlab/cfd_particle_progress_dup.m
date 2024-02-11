function cfd_particle_progress_dup(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 1.0; % work well for remote and local for dual screen

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

accruedTime = evalin('caller', 'accruedTime');
avgDen      = evalin('caller', 'avgDen');
avgVel      = evalin('caller', 'avgVel');
avgPrs      = evalin('caller', 'avgPrs');
avgVelGap   = evalin('caller', 'avgVelGap');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

subplot(1,4,1)
forceH = plot(accruedTime, avgDen, 'r-*', 'LineWidth', 2);
%title('average density inside particle volume');
xlabel('time (s)');
ylabel('inside avgDensity');
%legend('inside avgDensity');

subplot(1,4,2)
forceH = plot(accruedTime, avgVel, 'b-o', 'LineWidth', 2);
%title('average velocity inside particle volume');
xlabel('time (s)');
ylabel('inside avgVelocity');
%legend('inside avgVelocity');

subplot(1,4,3)
forceH = plot(accruedTime, avgPrs, 'm-x', 'LineWidth', 2);
%title('average pressure inside particle volume');
xlabel('time (s)');
ylabel('inside avgPressure');
%legend('inside avgPressure');

subplot(1,4,4)
forceH = plot(accruedTime, avgVelGap, 'g-d', 'LineWidth', 2);
%title('average pressure inside particle volume');
xlabel('time (s)');
ylabel('inside velocity gap');
%legend('inside velocity gap');

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '_dup.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_dup.png'), options);
