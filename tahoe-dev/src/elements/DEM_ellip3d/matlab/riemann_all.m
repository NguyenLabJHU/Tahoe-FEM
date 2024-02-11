fileInfo = dir('riemann_*.dat');
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1));
    riemann(fileInfo(i).name);
end
