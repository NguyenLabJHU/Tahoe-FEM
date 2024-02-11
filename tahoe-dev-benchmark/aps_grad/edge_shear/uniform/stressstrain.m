clear all
%
data1=load('eldata_mesh1_1.txt');
time1=data1(:,1);
gamma_x1=data1(:,2);
gamma_y1=data1(:,3);
effstr1=data1(:,4);
s_xz1=data1(:,5);
s_yz1=data1(:,6);
J2_1=1e-6*data1(:,7);
kappa1=1e-6*data1(:,8);
%
data10=load('eldata_mesh10_1.txt');
time10=data1(:,1);
gamma_x10=data1(:,2);
gamma_y10=data1(:,3);
effstr10=data10(:,4);
s_xz10=data10(:,5);
s_yz10=data10(:,6);
J2_10=1e-6*data10(:,7);
kappa10=1e-6*data10(:,8);
%
data2=load('out.dat');
time2=data2(:,1);
strain2=data2(:,3);
stress2=1e-6*data2(:,2);
kappa2=1e-6*data2(:,5);
%
data3=load('s_alpha.txt');
time3=data3(:,1);
strain3=1e-3*time3;
stress3=1e-6*data3(:,2);
%
data4=load('kappa.txt');
time4=data4(:,1);
strain4=1e-3*time4;
kappa4=1e-6*data4(:,2);
%
figure(1)
plot(strain2,stress2,gamma_x1,J2_1,gamma_x10,J2_10,strain3,stress3)
axis([0 0.3 0 300]);
xlabel('SHEAR STRAIN')
ylabel('SHEAR STRESS (MPa)')
legend('ODE solution','1 standard element','100 standard elements','100 vector elements')
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(2)
plot(strain2,kappa2,gamma_x1,kappa1,gamma_x10,kappa10,strain4,kappa4)
axis([0 0.3 0 90])
xlabel('SHEAR STRAIN')
ylabel('KAPPA (MPa)')
legend('ODE solution','1 standard element','100 standard elements','100 vector elements')
set(gca,'FontName','Helvetica','FontSize',16)
%
