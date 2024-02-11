close all;
clear all;

data1=load('three_d_hex27_8m_top_harmonic_run170_767.txt'); % top center
data2=load('three_d_hex27_8m_top_harmonic_run171_767.txt'); % top center 
data3=load('three_d_hex27_8m_top_harmonic_run172_767.txt'); % top center
data4=load('three_d_hex27_8m_top_harmonic_run170_769.txt'); % top corner
data5=load('three_d_hex27_8m_top_harmonic_run171_769.txt'); % top corner 
data6=load('three_d_hex27_8m_top_harmonic_run172_769.txt'); % top corner
data7=load('three_d_hex27_8m_top_harmonic_run170_768.txt'); % top corner
data8=load('three_d_hex27_8m_top_harmonic_run171_768.txt'); % top corner
data9=load('three_d_hex27_8m_top_harmonic_run172_768.txt'); % top corner
data10=load('three_d_hex27_8m_top_harmonic_run170_001.txt'); % bottom
data11=load('three_d_hex27_8m_top_harmonic_run171_001.txt'); % bottom
data12=load('three_d_hex27_8m_top_harmonic_run172_001.txt'); % bottom

time=data1(:,1);

d_y1=data1(:,7);
d_y2=data2(:,7);
d_y3=data3(:,7);
d_y4=data4(:,7);
d_y5=data5(:,7);
d_y6=data6(:,7);
d_y7=data7(:,7);
d_y8=data8(:,7);
d_y9=data9(:,7);
d_y10=data10(:,7);
d_y11=data11(:,7);
d_y12=data12(:,7);



p1=data1(:,15);
p2=data2(:,15);
p3=data3(:,15);
p4=data4(:,15);
p5=data5(:,15);
p6=data6(:,15);
p7=data7(:,15);
p8=data8(:,15);
p9=data9(:,15);
p10=data10(:,15);
p11=data11(:,15);
p12=data12(:,15);


s22_1=data1(:,10);
s22_2=data2(:,10);
s22_3=data3(:,10);
s22_4=data4(:,10);
s22_5=data5(:,10);
s22_6=data6(:,10);
s22_7=data7(:,10);
s22_8=data8(:,10);
s22_9=data9(:,10);
s22_10=data10(:,10);
s22_11=data11(:,10);
s22_12=data12(:,10);



ns_1=data1(:,22);
ns_2=data2(:,22);
ns_3=data3(:,22);
ns_4=data4(:,22);
ns_5=data5(:,22);
ns_6=data6(:,22);
ns_7=data7(:,22);
ns_8=data8(:,22);
ns_9=data9(:,22);
ns_10=data10(:,22);
ns_11=data11(:,22);
ns_12=data12(:,22);



nf_1=data1(:,23);
nf_2=data2(:,23);
nf_3=data3(:,23);
nf_4=data4(:,23);
nf_5=data5(:,23);
nf_6=data6(:,23);
nf_7=data7(:,23);
nf_8=data8(:,23);
nf_9=data9(:,23);
nf_10=data10(:,23);
nf_11=data11(:,23);
nf_12=data12(:,23);

%{
figure(1);
plot(time,d_y1,'b.-','MarkerSize',15)
hold on;
plot(time,d_y2,'r.-','MarkerSize',11)
hold on;
plot(time,d_y3,'g.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(meters)')
legend('k=0.1m/s','k=0.01m/s','k=0.001m/s')
text(0.05,-0.22,'NODE A');

figure(2);
plot(time,d_y4,'b.-','MarkerSize',15)
hold on;
plot(time,d_y5,'r.-','MarkerSize',11)
hold on;
plot(time,d_y6,'g.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(meters)')
legend('k=0.1m/s','k=0.01m/s','k=0.001m/s')
text(0.05,-0.15,'NODE B');

figure(3);
plot(time,d_y7,'b.-','MarkerSize',15)
hold on;
plot(time,d_y8,'r.-','MarkerSize',11)
hold on;
plot(time,d_y9,'g.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(meters)')
legend('k=0.1m/s','k=0.01m/s','k=0.001m/s')
text(0.05,-0.10,'NODE C');
%}

figure(4);
plot(time,d_y1,'b.-','MarkerSize',7)
hold on;
plot(time,d_y4,'b.-','MarkerSize',7)
hold on;
plot(time,d_y7,'b.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(meters)')
legend('k=0.1m/s')
text(0.05,-0.50,'NODE A');
text(0.05,-0.195,'NODE B');
text(0.05,-0.12,'NODE C');

figure(5);
plot(time,d_y2,'b.-','MarkerSize',7)
hold on;
plot(time,d_y5,'b.-','MarkerSize',7)
hold on;
plot(time,d_y8,'b.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(meters)')
legend('k=0.01m/s')
text(0.05,-0.45,'NODE A');
text(0.05,-0.18,'NODE B');
text(0.05,-0.12,'NODE C');

figure(6);
plot(time,d_y3,'b.-','MarkerSize',7)
hold on;
plot(time,d_y6,'b.-','MarkerSize',7)
hold on;
plot(time,d_y9,'b.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(meters)')
legend('k=0.001m/s')
text(0.05,-0.43,'NODE A');
text(0.05,-0.16,'NODE B');
text(0.05,-0.10,'NODE C');




figure(7);
plot(time,p10./1000,'b.-','MarkerSize',15)
hold on;
plot(time,p11./1000,'r.-','MarkerSize',11)
hold on;
plot(time,p12./1000,'g.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('EXCESS PORE PRESSURE(kPa)')
legend('k=0.1m/s','k=0.01m/s','k=0.001m/s')
text(0.01,-45,'NODE D')


figure(8);
plot(time,s22_10./1000,'b.-','MarkerSize',15)
hold on;
plot(time,s22_11./1000,'r.-','MarkerSize',11)
hold on;
plot(time,s22_12./1000,'g.-','MarkerSize',7)
hold on;
xlabel('TIME(seconds)')
ylabel('EFFECTIVE STRESS(kPa)')
legend('k=0.1m/s','k=0.01m/s','k=0.001m/s')
text(0.35,-10,'NODE D')


