function [d,u,m,p,e] = riemann_point(fileToRead1, position)

newData1 = importdata(fileToRead1);
numeric = newData1.data;

point = numeric( numeric(:, 1) == position, :);

colheaders = {'position' 'density' 'velocity' 'momentum' 'pressure' 'int_energy'};
for i = 1:size(point, 2)
    assignin('caller', genvarname(colheaders{i}), point(:,i));
end

position = evalin('caller', 'position');
density  = evalin('caller', 'density');
velocity = evalin('caller', 'velocity');
momentum = evalin('caller', 'momentum');
pressure = evalin('caller', 'pressure');
int_energy = evalin('caller', 'int_energy');

d = density;
u = velocity;
m = momentum;
p = pressure;
e = int_energy;
