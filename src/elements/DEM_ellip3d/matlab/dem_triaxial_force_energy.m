function dem_triaxial_force_energy(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 1.0; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

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

subplot(1,2,1)
[AX,H1,H2] = plotyy(epsilon_z, [avgNormal, avgShear], epsilon_z, avgPenetr); % such that avgNormal and avgShear show all y range
grid on;
grid minor;
pbaspect(AX(1), [1 1 1]);
pbaspect(AX(2), [1 1 1]);
%set(H1, 'Color', 'b', 'LineWidth', 3);
set(H1, 'LineWidth', 3);
set(H2, 'Color', 'r', 'LineWidth', 3);
set(AX(1),'YColor','b');
set(AX(2),'YColor','r');
%line(epsilon_z, avgShear, 'Parent', AX(1), 'Color', 'k', 'LineWidth', 3); % unnecessary
title('Force and penetration');
xlabel('Axial strain'); 
set(get(AX(1), 'Ylabel'), 'String', 'Force (N)');
set(get(AX(2), 'Ylabel'), 'String', 'Penetration (m)');
%set(AX(1), 'xlim', [0, 3.5e-4]);
%set(AX(2), 'xlim', [0, 4000]);
lh=legend('avgNormal', 'avgShear', 'avgPenetr', 'location', 'east');
set(lh, 'box', 'off');

subplot(1,2,2)
%forceH = plot(accruedTime, normal_z1, 'r--o', accruedTime, normal_z2, 'm--', accruedTime, normal_x1, 'b-x', accruedTime, normal_x2, 'k-p','LineWidth', 1); 
forceH = plot(epsilon_z, transEnergy, epsilon_z, rotatEnergy, epsilon_z, kinetEnergy, 'LineWidth', 3); 
grid on;
grid minor;
pbaspect([1 1 1])
title('Energy');
xlabel('Axial strain');
ylabel('Energy (J)');
lh=legend('transEnergy', 'rotatEnergy', 'kinetEnergy', 'location', 'east');
set(lh, 'box', 'off');
%axis([0, 3.5e-4, ylim]);
%xlim([0, 4000]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_force.png'), options);
