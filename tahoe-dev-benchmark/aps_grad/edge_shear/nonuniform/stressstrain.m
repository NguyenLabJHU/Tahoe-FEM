clear all
%
data1=load('R_curl_10.txt');
time1=data1(:,1);
displ1=time1*1e-3;
react1=1e-3*data1(:,2);
%
data2=load('R_nocurl_10.txt');
time2=data2(:,1);
displ2=time2*1e-3;
react2=1e-3*data2(:,2);
%
data3=load('R_nocurl_vect.txt');
displ3=1e3*data3(:,1);
react3=1e-3*data3(:,2);
%
figure(1)
plot(displ1,react1,displ2,react2,displ3,react3)
axis([0 0.3 0 250 ]);
xlabel('DISPLACEMENT (mm)')
%xlabel('time')
ylabel('REACTION (kN)')
legend('100 standard elements: L=0.1mm','100 standard elements: L=0', '100 vector elements: L=0')
set(gca,'FontName','Helvetica','FontSize',16)
%


