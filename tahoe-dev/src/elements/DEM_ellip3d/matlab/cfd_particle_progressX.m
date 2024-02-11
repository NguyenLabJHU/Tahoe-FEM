function cfd_particle_progress(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

accruedTime = evalin('caller', 'accruedTime');
penalFx     = evalin('caller', 'penalFx');
pressureFx  = evalin('caller', 'pressureFx');
internalFx  = evalin('caller', 'internalFx');
totalFx     = penalFx + pressureFx + internalFx;
accelX      = evalin('caller', 'accelX');
velocX      = evalin('caller', 'velocX');
accelY      = evalin('caller', 'accelY');
velocY      = evalin('caller', 'velocY');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'off');

subplot(1,2,1)
forceH = plot(accruedTime, penalFx, 'r--o', accruedTime, pressureFx, 'm--', accruedTime, internalFx, 'b-x', accruedTime, totalFx, 'k-p','LineWidth', 2); 
title('drag force');
xlabel('time (s)');
ylabel('force (N)');
legend('viscousFx', 'pressureFx', 'internalFx', 'totalFx'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);

subplot(1,2,2)
[AX,H1,H2] = plotyy(accruedTime, accelX, accruedTime, velocX);
set(H1, 'Color', 'b', 'LineStyle', '--', 'Marker', 'o', 'LineWidth', 2);
set(H2, 'Color', 'r', 'LineStyle', '-',  'Marker', 'x', 'LineWidth', 2);
set(AX(1),'YColor','b');
set(AX(2),'YColor','r');
title('kinematics');
xlabel('time (s)'); 
set(get(AX(1), 'Ylabel'), 'String', 'acceleration (m/s^2)');
set(get(AX(2), 'Ylabel'), 'String', 'velocity (m/s)');
%set(AX(1), 'xlim', [0, 3.5e-4]);
%set(AX(2), 'xlim', [0, 0.4e-4]);
legend('accelX', 'velocX'); % , 'location', 'best'

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(strcat(fileToRead1, '_X'), '.png'), options);
