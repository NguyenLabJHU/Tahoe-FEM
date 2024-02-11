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

[AX,H1,H2] = plotyy(dia, number, dia, mass);
%pbaspect([2 1 1]);
%set(AX,'PlotBoxAspectRatio',[1.5 1 1]);
set(H1, 'Color', 'b', 'LineStyle', '--', 'Marker', 'o', 'MarkerSize', 12, 'LineWidth', 4);
set(H2, 'Color', 'r', 'LineStyle', '-',  'Marker', 'x', 'MarkerSize', 12, 'LineWidth', 4);
set(AX(1),'YColor','b');
set(AX(2),'YColor','r');
set(AX(1),'YTick',0:10:100);
set(AX(2),'YTick',0:10:100);

set(AX,'XScale','log');
custom_xticks = fliplr(dia');
set(AX,'XTick',custom_xticks);
xtickangle(60);
title('Number vs mass percentage');
xlabel('Particle size (mm, diameter)'); 
set(get(AX(1), 'Ylabel'), 'String', 'Number percentage smaller');
set(get(AX(2), 'Ylabel'), 'String', 'Mass percentage smaller');
%set(AX(1), 'xlim', [0, 3.5e-4]);
%set(AX(2), 'xlim', [0, 0.4e-4]);
legend('number percentage', 'mass percentage', 'location', 'NorthWest');

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(AX, 'fontsize', 18);
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
options.Format = 'png';
[pathstr,name,ext] = fileparts(fileToRead1);
hgexport(fh, strcat(name, '.png'), options);
