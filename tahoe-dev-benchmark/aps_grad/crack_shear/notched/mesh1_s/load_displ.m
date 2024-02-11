clear all
%
data1=load('R_nocurl.txt');
displ1=1e3*abs(data1(:,1));
R1=1e-3*data1(:,2);
%
data2=load('R_curl.txt');
displ2=1e3*abs(data2(:,1));
R2=1e-3*data2(:,2);
%
data3=load('R_nocurl_zero.txt');
displ3=1e3*abs(data3(:,1));
R3=1e-3*data3(:,2);
%
data4=load('R_curl_zero.txt');
displ4=1e3*abs(data4(:,1));
R4=1e-3*data4(:,2);
%
%
figure(1)
plot(displ4,R4,displ3,R3,displ2,R2,displ1,R1,'LineWidth',1)
xlabel('DISPLACEMENT (mm)')
ylabel('REACTION (kN)')
legend('L = 0.1mm, \gamma^p = 0 on boundary','L = 0, \gamma^p = 0 on boundary','L = 0.1 mm','L = 0')
set(gca,'FontName','Helvetica','FontSize',16)
%
figure(2)
plot(displ2,R2,displ1,R1,'LineWidth',1)
xlabel('DISPLACEMENT (mm)')
ylabel('REACTION (kN)')
legend('L = 0.1 mm','L = 0')
set(gca,'FontName','Helvetica','FontSize',16)
%
