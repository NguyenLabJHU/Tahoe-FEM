function tuneMass(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % 0.5 works for remote and local dual screens; 1.0 works for remote and local single screen.

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

number = 100.0 * evalin('caller', 'number_percent');
mass   = 100.0 * evalin('caller', 'mass_percent');
dia    = 1.0E+3 * evalin('caller', 'diameter');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'on');

plot(dia, number, 'b--o', dia, mass, 'r-x', 'MarkerSize', 12, 'LineWidth', 4);
%pbaspect([2 1 1]);

yticks(0:10:100);
set(gca,'XScale','log');
custom_xticks = fliplr(dia');
set(gca,'XTick',custom_xticks);
xtickangle(60);
title('Number vs mass percentage');
xlabel('Particle size (mm, diameter)'); 
ylabel('Percentage smaller');
legend('number percentage', 'mass percentage', 'location', 'SouthEast'); %NorthWest
xlim([0.2, 2]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
%set(gca, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
[pathstr,name,ext] = fileparts(fileToRead1);
hgexport(fh, strcat(name, '.png'), options);
