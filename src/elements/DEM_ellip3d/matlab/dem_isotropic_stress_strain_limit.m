function dem_isotropic_stress(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

timeStep = 1; %5e-7;
newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

step = evalin('caller', 'iteration');
traction_x1 = evalin('caller', 'traction_x1');
traction_x2 = evalin('caller', 'traction_x2');
traction_y1 = evalin('caller', 'traction_y1');
traction_y2 = evalin('caller', 'traction_y2');
traction_z1 = evalin('caller', 'traction_z1');
traction_z2 = evalin('caller', 'traction_z2');
mean_stress  = evalin('caller', 'mean_stress ');
bulk_volume = evalin('caller', 'bulk_volume');
density = evalin('caller', 'density');
epsilon_x = evalin('caller', 'epsilon_x');
epsilon_y = evalin('caller', 'epsilon_y');
epsilon_z = evalin('caller', 'epsilon_z');
epsilon_v = evalin('caller', 'epsilon_v');
void_ratio = evalin('caller', 'void_ratio');
porosity = evalin('caller', 'porosity');
avgNormal = evalin('caller', 'avgNormal');
avgShear = evalin('caller', 'avgShear');
avgPenetr = evalin('caller', 'avgPenetr');
transEnergy = evalin('caller', 'transEnergy');
rotatEnergy = evalin('caller', 'rotatEnergy');
kinetEnergy = evalin('caller', 'kinetEnergy');
accruedTime = step * timeStep;
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'on');

%{
subplot(1,2,1)
%forceH = plot(accruedTime, traction_z1, 'r--o', accruedTime, traction_z2, 'm--', accruedTime, traction_x1, 'b-x', accruedTime, traction_x2, 'k-p','LineWidth', 1); 
forceH = plot(accruedTime, traction_z2, 'r', accruedTime, traction_z1, 'b', accruedTime, traction_x1, 'g', accruedTime, traction_x2, 'k', accruedTime, traction_y1, 'c', accruedTime, traction_y2, 'm','LineWidth', 3); 
grid on;
grid minor;
pbaspect([1 1 1]);
title('Traction on 6 walls');
xlabel('Step');
ylabel('Traction (N/m^2)');
lh=legend('top', 'bottom', 'front', 'back', 'left', 'right', 'location', 'northwest');
set(lh, 'box', 'off');
%axis([0, 3.5e-4, ylim]);
ylim([1E+5, 7E+5]);
xlim([0, 4E+4]);
%}

subplot(1,1,1)
%forceH = plot(accruedTime, traction_z1, 'r--o', accruedTime, traction_z2, 'm--', accruedTime, traction_x1, 'b-x', accruedTime, traction_x2, 'k-p','LineWidth', 1); 
forceH = plot(epsilon_v, mean_stress, 'b', 'LineWidth', 3); 
grid on;
grid minor;
pbaspect([1 1 1]);
title('Stress vs strain');
xlabel('Volumetric strain');
ylabel('Mean stress (N/m^2)');
%legend('stress-strain', 'location', 'northwest');
%axis([0, 3.5e-4, ylim]);
ylim([1E+5, 7E+5]);
xlim([0, 0.01]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_stress.png'), options);
