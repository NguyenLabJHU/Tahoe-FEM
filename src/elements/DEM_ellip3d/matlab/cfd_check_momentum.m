%userpath('/home/yanb/matlab/')
screenWidth = 0.5; % work well for remote and local for dual screen

newData1 = importdata('momentum_progress');
% Create new variables in the caller workspace from those fields.
for i = 1:size(newData1.colheaders, 2)
    assignin('base', genvarname(newData1.colheaders{i}), newData1.data(:,i));
end

fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);

subplot(1,2,1)
plot(iteration, momAFluidZ, 'r--o', iteration, momBFluidZ, 'b-x', iteration, momPtclZ, 'g:s', ...
     iteration, momAZ, 'c-*', iteration, momBZ, 'k-d', 'LineWidth', 2); % o x
title('momentum');
xlabel('iteration');
ylabel('momentum');
legend('momAFluidZ', 'momBFluidZ', 'momPtclZ', 'momAZ', 'momBZ', 'location', 'best'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);

subplot(1,2,2)
plot(iteration, momPtclZ, 'g:s', iteration, momAZ, 'c-*', iteration, momBZ, 'k-d', 'LineWidth', 2); % o x
title('momentum');
xlabel('iteration');
ylabel('momentum');
legend('momPtclZ', 'momAZ', 'momBZ'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);

%subplot(1,3,3)
%plot(iteration, momDiffZ, 'm--o', 'LineWidth', 2); % o x
%title('momentum');
%xlabel('iteration');
%ylabel('momentum');
%legend('momDiffZ'); %, 'location', 'best'
%axis([0, 3.5e-4, ylim]);
%xlim([0, 3.5e-4]);

set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, 'momentum_progress.png', 'png');
options.Format = 'png';
hgexport(fh, 'momentum_progress.png', options);
