clear all
%
data1=load('R.txt');
time1=data1(:,1);
time1=time1*(1e-3);
R1=data1(:,2);
%
data2=load('R_curl.txt');
time2=data2(:,1);
time2=time2*(1e-3);
R2=data2(:,2);
%
data3=load('R_zero.txt');
time3=data3(:,1);
time3=time3*(1e-3);
R3=data3(:,2);
%
data4=load('R_curl_zero.txt');
time4=data4(:,1);
time4=time4*(1e-3);
R4=data4(:,2);
%
figure(1)
plot(time1,R1,time2,R2)
xlabel('DISPLACEMENT (mm)')
ylabel('REACTION (N)')
legend('without curl','with curl')
%
figure(2)
plot(time1,R1,time2,R2,time3,R3,time4,R4)
xlabel('DISPLACEMENT (mm)')
ylabel('REACTION (N)')
legend('without curl','with curl','without curl: gammap=0','with curl: gammap=0')
%



