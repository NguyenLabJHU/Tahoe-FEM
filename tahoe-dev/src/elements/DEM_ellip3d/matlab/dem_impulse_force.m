function dem_impulse_force(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

timeStep = 1; %5e-7;
newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

step = evalin('caller', 'iteration');
normal_x1 = evalin('caller', 'normal_x1');
normal_x2 = evalin('caller', 'normal_x2');
normal_y1 = evalin('caller', 'normal_y1');
normal_y2 = evalin('caller', 'normal_y2');
normal_z1 = evalin('caller', 'normal_z1');
normal_z2 = evalin('caller', 'normal_z2');
transEnergy = evalin('caller', 'transEnergy');
rotatEnergy = evalin('caller', 'rotatEnergy');
accruedTime = step * timeStep;
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'on');

subplot(1,1,1)
%forceH = plot(accruedTime, normal_z1, 'r--o', accruedTime, normal_z2, 'm--', accruedTime, normal_x1, 'b-x', accruedTime, normal_x2, 'k-p','LineWidth', 1); 
forceH = plot(accruedTime, normal_z2, 'r', accruedTime, normal_z1, 'b', accruedTime, normal_x1, 'g', accruedTime, normal_x2, 'k', accruedTime, normal_y1, 'c', accruedTime, normal_y2, 'm','LineWidth', 3); 
grid on;
grid minor;
pbaspect([1 1 1]);
title('Wall response to a top impulse');
xlabel('Step');
ylabel('Force (N)');
lh=legend('top', 'bottom', 'front', 'back', 'left', 'right'); %, 'location', 'best');
set(lh, 'box', 'off');
%axis([0, 3.5e-4, ylim]);
xlim([0, 4000]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_force.png'), options);
