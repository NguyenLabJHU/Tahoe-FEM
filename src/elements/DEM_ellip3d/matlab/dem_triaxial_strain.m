function dem_triaxial_strain(fileToRead1)
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
sigma1 = (traction_z1 + traction_z2)/2.0;
sigma3 = (traction_x1 + traction_x2 + traction_y1 + traction_y2)/4.0;

Green_x = evalin('caller', 'Green_x');
Green_y = evalin('caller', 'Green_y');
Green_z = evalin('caller', 'Green_z');
Euler_x = evalin('caller', 'Euler_x');
Euler_y = evalin('caller', 'Euler_y');
Euler_z = evalin('caller', 'Euler_z');
vol_strain = evalin('caller', 'vol_strain');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'on');

subplot(1,2,1)
%forceH = plot(accruedTime, traction_z1, 'r--o', accruedTime, traction_z2, 'm--', accruedTime, traction_x1, 'b-x', accruedTime, traction_x2, 'k-p','LineWidth', 1); 
forceH = plot(step, -epsilon_z, 'r', step, Green_z, 'b', step, Euler_z, 'g', step, -epsilon_x, 'k', step, Green_x, 'c', step, Euler_x, 'm','LineWidth', 3); 
grid on;
grid minor;
pbaspect([1 1 1]);
title('Various strain measures');
xlabel('Steps');
ylabel('Strain');
lh=legend('small_z', 'Green_z', 'Euler_z', 'small_x', 'Green_x', 'Euler_x', 'location', 'southwest');
set(lh, 'box', 'off');
%axis([0, 3.5e-4, ylim]);
%xlim([0, 4000]);

subplot(1,2,2)
%forceH = plot(accruedTime, traction_z1, 'r--o', accruedTime, traction_z2, 'm--', accruedTime, traction_x1, 'b-x', accruedTime, traction_x2, 'k-p','LineWidth', 1); 
forceH = plot(step, -epsilon_v, 'r', step, vol_strain, 'b', 'LineWidth', 3); 
grid on;
grid minor;
pbaspect([1 1 1]);
title('Volume strain measures');
xlabel('Steps');
ylabel('Volume starin');
legend('small strain', 'finite strain', 'location', 'south');
%axis([0, 3.5e-4, ylim]);
%xlim([0, 4000]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_strain.png'), options);
