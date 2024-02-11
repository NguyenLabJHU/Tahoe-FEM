clear all
%
% node 1
sim=load('ndata_1.txt');
time=(sim(:,1));
d=(sim(:,2));
sig33=(sim(:,3));
c_f_ref=(sim(:,4));
c_f_cur=(sim(:,5));
%
sim=load('ndata_1_x10e5.txt');
time_a=(sim(:,1));
d_a=(sim(:,2));
sig33_a=(sim(:,3));
c_f_ref_a=(sim(:,4));
c_f_cur_a=(sim(:,5));
%
figure(1)
plot(time,abs(sig33), time_a, abs(sig33_a))
xlabel('time (sec)')
ylabel('stress (Pa)')
legend('D = 3.06e-12 sec', 'D = 3.06e-5 sec') 
title('node 1')
%
figure(2)
plot(time,c_f_ref, time_a, c_f_ref_a)
xlabel('time (sec)')
ylabel('ref concentration (kg/m^3)')
legend('D = 3.06e-12 sec', 'D = 3.06e-5 sec') 
title('node 1')
%
