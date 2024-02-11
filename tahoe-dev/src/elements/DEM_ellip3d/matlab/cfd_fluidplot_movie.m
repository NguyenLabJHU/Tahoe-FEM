fileInfo = dir('couple_fluidplot_*.png');
writerObj = VideoWriter('couple_fluidplot.avi', 'Motion JPEG AVI');
writerObj.FrameRate = 4;
open(writerObj);
for i = 1 : size(fileInfo, 1)
  filename = sprintf('couple_fluidplot_%03d.dat.png', i-1);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);
