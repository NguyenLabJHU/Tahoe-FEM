function dem_contact(fileToRead1)

% 1-scattered
% 2-surface mesh
% 3-hull
type = 3;
zoom = 0;
scatterSize = 0; %only work for zoom==1: 0-equal size; 1-proportional size

newData1 = importdata(fileToRead1, ' ', 3); %delimiterIn,headerlinesIn
%disp(newData1);
numeric  = newData1.data;

x = numeric(:, 1);
y = numeric(:, 2);
z = numeric(:, 3);
normal = numeric(:, 13);
shear  = numeric(:, 14);
total  = numeric(:, 15);

screenWidth = 1.0;
fh = figure('units', 'normalized', 'outerposition', [0 0 screenWidth 1]);
set(fh, 'visible', 'on');

if type == 1
 if scatterSize == 1
  scatter3(x,y,z,total*1000,'r','filled');
 elseif scatterSize == 0
  scatter3(x,y,z,10,'r','filled');
 end
elseif type == 2
  tri = delaunay(x,y,z);
  trimesh(tri,x,y,z);
  size(tri)
elseif type == 3
  dt = delaunayTriangulation(x,y,z);
  [K, volume]=convexHull(dt);
  volume
  %trisurf(K,dt.Points(:,1),dt.Points(:,2),dt.Points(:,3))
  tetramesh(dt); % meshes in type 2 and type 3 are the same?
  %camorbit(20,0)
end

if zoom == 0
  xlim([0, 0.5]);
  ylim([0, 0.5]);
  zlim([0, 0.2]);
elseif zoom == 1
  xlim([0.2, 0.3]);
  ylim([0.2, 0.3]);
  zlim([0, 0.12]);
end

set(gca,'DataAspectRatio',[1 1 1]);
set(findall(gcf, '-property', 'fontSize'), 'fontSize', 26, 'fontWeight', 'bold');
set(gcf, 'paperpositionmode', 'auto');
%saveas(fh, strcat(fileToRead1, '.png'), 'png');
%options.Format = 'png';
%hgexport(fh, strcat(fileToRead1, '_force.png'), options);
%print(strcat(fileToRead1, '.png'), '-dpng', '-r600');

%view(90, 90);
%export_fig(strcat(fileToRead1, '_top.png'), '-r600');
%view(90, -90);
%export_fig(strcat(fileToRead1, '_bottom.png'), '-r600');

view(90, 0);
%view(3);

if zoom == 0
 if type == 1
   fout=strcat(fileToRead1, '_scatter.png');
 elseif type == 2
   fout=strcat(fileToRead1, '_mesh.png');
 elseif type == 3
   fout=strcat(fileToRead1, '_hull.png');
 end
elseif zoom == 1
 if type == 1
  if scatterSize == 1
   fout=strcat(fileToRead1, '_scatter_zoom.png');
  elseif scatterSize == 0
   fout=strcat(fileToRead1, '_scatter_zoom_equal.png');
  end
 elseif type == 2
   fout=strcat(fileToRead1, '_mesh_zoom.png');
 elseif type == 3
   fout=strcat(fileToRead1, '_hull_zoom.png');
 end
end
%export_fig(fout, '-r600');


