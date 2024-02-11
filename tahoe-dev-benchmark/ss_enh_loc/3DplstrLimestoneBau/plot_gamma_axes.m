function dummy = plot_gamma_axes(center, radius, psi, handle);
dummy = 1;

%center = [0];
%radius = 1;
%psi = 0.8;
%handle = 1;

%Plot_Gamma(center, radius, 1, handle);
%hold on;
Plot_Gamma(center, radius, psi, handle);
hold on;
plot([0;0],[0;1.*radius],'k-');
plot([0; .625*sqrt(3)*radius], [0;-.625*radius],'k-');
plot([0; -.625*sqrt(3)*radius], [0;-.625*radius],'k-');