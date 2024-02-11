function dem_deposit_force_energy(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

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
kinetEnergy = evalin('caller', 'kinetEnergy');
accruedTime = step;
normal_side = (normal_x1 + normal_x2 + normal_y1 +normal_y2) / 4.0;
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'on');

subplot(1,1,1)
[AX,H1,H2] = plotyy(accruedTime, [normal_z1, normal_side], accruedTime, [transEnergy, rotatEnergy, kinetEnergy]);
grid on;
grid minor;
pbaspect(AX(1), [1 1 1]);
pbaspect(AX(2), [1 1 1]);
set(H1, 'LineWidth', 3);
set(H2, 'LineWidth', 3);
set(AX(1),'YColor','b');
set(AX(2),'YColor','r');
%line(accruedTime, normal_z1, 'Parent', AX(1));
%line(accruedTime, rotatEnergy, 'Parent', AX(2));
title('Force and energy');
xlabel('Step'); 
set(get(AX(1), 'Ylabel'), 'String', 'Force (N)');
set(get(AX(2), 'Ylabel'), 'String', 'Energy (J)');
%set(AX(1), 'xlim', [0, 3.5e-4]);
%set(AX(2), 'xlim', [0, 4000]);
lh=legend('bottomForce', 'sideForce', 'transEnergy', 'rotatEnergy', 'kinetEnergy'); % , 'location', 'best'
set(lh, 'box', 'off');

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_force_energy.png'), options);
