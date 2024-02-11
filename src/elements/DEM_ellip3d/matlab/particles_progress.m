function particles_progress(fileToRead1)
screenWidth = 0.5; % work well for remote and local for dual screen

newData1 = importdata(fileToRead1);

for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

accruedTime = evalin('caller', 'accruedTime');
avgNormal = evalin('caller', 'avgNormal');
avgShear = evalin('caller', 'avgShear');
transEnergy = evalin('caller', 'transEnergy');
rotatEnergy = evalin('caller', 'rotatEnergy');
normal_x1 = evalin('caller', 'normal_x1');
normal_x2 = evalin('caller', 'normal_x2');
normal_y1 = evalin('caller', 'normal_y1');
normal_y2 = evalin('caller', 'normal_y2');
normal_z1 = evalin('caller', 'normal_z1');
normal_z2 = evalin('caller', 'normal_z2');
contact_x1 = evalin('caller', 'contact_x1');
contact_x2 = evalin('caller', 'contact_x2');
contact_y1 = evalin('caller', 'contact_y1');
contact_y2 = evalin('caller', 'contact_y2');
contact_z1 = evalin('caller', 'contact_z1');
contact_z2 = evalin('caller', 'contact_z2');
contact_inside = evalin('caller', 'contact_inside');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1.0]);

subplot(1,4,1)
forceH = plot(accruedTime, avgNormal, 'r-', accruedTime, avgShear, 'b-', 'LineWidth', 2); % o x
title('internal force');
xlabel('time');
ylabel('force');
legend('normal', 'shear', 'location', 'best'); %, 'location', 'best'
axis([0, 0.5, ylim]);

subplot(1,4,2)
forceH = plot(accruedTime, (contact_x1+contact_x2)/2, 'r-', accruedTime, (contact_y1+contact_y2)/2, 'b-', accruedTime, contact_z1, 'k-', accruedTime, contact_inside, 'g-', 'LineWidth', 2);
title('contact');
xlabel('time');
ylabel('contact');
legend('x', 'y', 'z1', 'inside', 'location', 'best');
axis([0, 0.5, ylim]);

subplot(1,4,3)
energyH = plot(accruedTime, transEnergy, 'm-', accruedTime, rotatEnergy, 'g-', 'LineWidth', 2);
title('energy');
xlabel('time');
ylabel('energy');
legend('transEnergy', 'rotateEnergy', 'location', 'best');
axis([0, 0.5, ylim]);

subplot(1,4,4)
forceH = plot(accruedTime, (normal_x1+normal_x2)/2, 'r-', accruedTime, (normal_y1+normal_y2)/2, 'b-', accruedTime, normal_z1, 'k-', 'LineWidth', 2);
title('boundary force');
xlabel('time');
ylabel('force');
legend('x', 'y', 'z1', 'location', 'best');
axis([0, 0.5, ylim]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 18, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
saveas(fh, 'particles_progress.png', 'png');
%options.Format = 'png';
%hgexport(fh, 'particles_progress.png', options);
