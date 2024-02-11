function cfd_particle_progress_Cd(fileToRead1)
%userpath('/home/yanb/matlab/')
screenWidth = 1.0; % work well for remote and local for dual screen

newData1 = importdata(fileToRead1);
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('caller', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

accruedTime = evalin('caller', 'accruedTime');
penalFz     = evalin('caller', 'penalFz');
pressureFz  = evalin('caller', 'pressureFz');
internalFz  = evalin('caller', 'internalFz');
totalFz     = penalFz + pressureFz + internalFz;
viscousCd   = evalin('caller', 'viscousCd');
pressureCd  = evalin('caller', 'pressureCd');
internalCd  = evalin('caller', 'internalCd');
totalCd     = evalin('caller', 'totalCd');

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

subplot(1,2,1)
forceH = plot(accruedTime, penalFz, 'r--o', accruedTime, pressureFz, 'm--', accruedTime, internalFz, 'b-x', accruedTime, totalFz, 'k-p','LineWidth', 2); 
title('drag force');
xlabel('time (s)');
ylabel('force (N)');
legend('viscousFz', 'pressureFz', 'internalFz', 'totalFz'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);

subplot(1,2,2)
CdH = plot(accruedTime, viscousCd, 'r--o', accruedTime, pressureCd, 'm--', accruedTime, internalCd, 'b-x', accruedTime, totalCd, 'k-p', 'LineWidth', 2); % o x
title('drag coefficient');
xlabel('time (s)');
ylabel('Cd');
legend('viscousCd', 'pressureCd', 'internalCd', 'totalCd'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);
ylim([0, 1.5]);
yticks(0:0.1:1.5);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '_Cd.png'), 'png');
options.Format = 'png';
hgexport(fh, strcat(fileToRead1, '_Cd.png'), options);
