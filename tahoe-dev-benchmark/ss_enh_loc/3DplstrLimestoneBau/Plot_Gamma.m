%function Plot_Gamma plots a circle on the complex plane with 
% center at center and radius or radius. Plot on the figure handle
function Plot_Gamma(center, radius, psi, handle);

xcoords = zeros(201);
ycoords = zeros(201);

for i = 0:200
   beta = i*pi/100;
   Gamma = .5*( 1 + sin(3*beta) + (1.0/psi) * (1 - sin(3*beta)));
   xcoords(i+1) = real(center) + (1/Gamma)*radius*cos(beta);
   ycoords(i+1) = imag(center) + (1/Gamma)*radius*sin(beta);
end

figure(handle);
plot(xcoords, ycoords,'b-');
%hold on;
%plot(xcoords, -1*ycoords);
