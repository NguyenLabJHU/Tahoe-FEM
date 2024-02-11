%program 1

x=0:0.1:2*pi;

[x,y]=meshgrid(x);

z=sin(y).*cos(x);

mesh(x,y,z);

xlabel('x-axis'),ylabel('y-axis'),zlabel('z-axis');

title('mesh'); pause;

%program 2

x=0:0.1:2*pi;

[x,y]=meshgrid(x);

z=sin(y).*cos(x);

surf(x,y,z);

xlabel('x-axis'),ylabel('y-axis'),zlabel('z-axis');

title('surf'); pause;

%program 3

x=0:0.1:2*pi;

[x,y]=meshgrid(x);

z=sin(y).*cos(x);

plot3(x,y,z);

xlabel('x-axis'),ylabel('y-axis'),zlabel('z-axis');

title('plot3-1');grid;