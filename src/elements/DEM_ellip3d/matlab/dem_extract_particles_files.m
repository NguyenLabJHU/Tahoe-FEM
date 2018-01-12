% 1. specify files
fileInfo = dir('deposit_particle_*');
for i = 1:size(fileInfo, 1)
    fprintf(1, 'step%4i of%4i\n', i, size(fileInfo, 1))
    dem_extract_particles(fileInfo(i).name);
end

% 2. specify ranges
% to change extraction range, modify dem_extract_particles.m
