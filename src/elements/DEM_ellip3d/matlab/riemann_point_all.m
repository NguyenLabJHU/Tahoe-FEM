timeSnap = 0.01;
position = 0.351;
suffix = '_pos0.351';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileInfo = dir('riemann_*.dat');
curve=zeros(size(fileInfo, 1), 6);
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    curve(i, 1) = i * timeSnap;
    [curve(i, 2), curve(i, 3), curve(i, 4), curve(i, 5), curve(i, 6)] = riemann_point(fileInfo(i).name, position);
end
%curve
%dlmwrite('point.dat', curve, 'delimiter', '\t', 'precision', 6);
fOut=fopen(strcat(strcat('point', suffix), '.dat'), 'w');
fprintf(fOut, '%15s%15s%15s%15s\n', 'time', 'position', 'density', 'velocity', 'momentum', 'pressure', 'int_energy');
for i = 1:size(curve,1)
    for j = 1:size(curve,2)
    fprintf(fOut, '%15.6e', curve(i, j));
    end
    fprintf(fOut, '\n');
end
fclose(fOut);

screenWidth = 0.5; % work well for remote and local for dual screen
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1.0]);

subplot(1,5,1)
plot(curve(:, 1), curve(:, 2), 'r-o', 'LineWidth', 2);
title('density');
xlabel('time');

subplot(1,5,2)
plot(curve(:, 1), curve(:, 3), 'r-o', 'LineWidth', 2);
title('velocity');
xlabel('time');

subplot(1,5,3)
plot(curve(:, 1), curve(:, 4), 'r-o', 'LineWidth', 2);
title('momentum');
xlabel('time');

subplot(1,5,4)
plot(curve(:, 1), curve(:, 5), 'r-o', 'LineWidth', 2);
title('pressure');
xlabel('time');

subplot(1,5,5)
plot(curve(:, 1), curve(:, 6), 'r-o', 'LineWidth', 2);
title('internal energy');
xlabel('time');

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
saveas(fh, strcat(strcat('point', suffix), '.png'), 'png');



