
%close all
clear all

format short e

sim=load('impulse_100kN_1em8_2577.txt');
time_100kN_1em8_2577=sim(:,1);
displ_100kN_1em8_2577=sim(:,4)*1000;

sim=load('impulse_100kN_2577.txt');
time_100kN_2577=sim(:,1);
displ_100kN_2577=sim(:,4)*1000;

% node displ time
figure(1)
plot(time_100kN_2577,displ_100kN_2577,'.-k',time_100kN_1em8_2577,displ_100kN_1em8_2577,'*-k','LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (mm)')
legend('solid','k^{hat} = 10^{-8} m^2/Pa.s')
axis([0 1 -1.5 1])
text(0.85,0.5,'NODE A','FontSize',16)
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('impulse_100kN_1em8_2457.txt');
time_100kN_1em8_2457=sim(:,1);
displ_100kN_1em8_2457=sim(:,4)*1000;

sim=load('impulse_100kN_2457.txt');
time_100kN_2457=sim(:,1);
displ_100kN_2457=sim(:,4)*1000;

% node displ time
figure(2)
plot(time_100kN_2457,displ_100kN_2457,'.-k',time_100kN_1em8_2457,displ_100kN_1em8_2457,'*-k','LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (mm)')
legend('solid','k^{hat} = 10^{-8} m^2/Pa.s')
axis([0 1 -0.8 0.6])
text(0.85,0.3,'NODE B','FontSize',16)
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})


sim=load('impulse_4MN_1em8_2577.txt');
time_4MN_1em8_2577=sim(:,1);
displ_4MN_1em8_2577=sim(:,4)*1000;

sim=load('impulse_4MN_2577.txt');
time_4MN_2577=sim(:,1);
displ_4MN_2577=sim(:,4)*1000;

% node displ time
figure(3)
plot(time_4MN_2577,displ_4MN_2577,'.-k',time_4MN_1em8_2577,displ_4MN_1em8_2577,'*-k','LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (mm)')
legend('solid','k^{hat} = 10^{-8} m^2/Pa.s')
axis([0 1 -80 60])
text(0.85,50,'NODE A','FontSize',16)
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('impulse_4MN_1em8_2457.txt');
time_4MN_1em8_2457=sim(:,1);
displ_4MN_1em8_2457=sim(:,4)*1000;

sim=load('impulse_4MN_2457.txt');
time_4MN_2457=sim(:,1);
displ_4MN_2457=sim(:,4)*1000;

% node displ time
figure(4)
plot(time_4MN_2457,displ_4MN_2457,'.-k',time_4MN_1em8_2457,displ_4MN_1em8_2457,'*-k','LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (mm)')
legend('solid','k^{hat} = 10^{-8} m^2/Pa.s')
axis([0 1 -40 30])
text(0.85,15,'NODE B','FontSize',16)
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('impulse_100kN_1em8_2577.txt');
time_100kN_1em8_2577=sim(:,1);
displ_100kN_1em8_2577=sim(:,4)*1000*40;

% node displ time
figure(5)
plot(time_100kN_1em8_2577,displ_100kN_1em8_2577,'.-k',time_4MN_1em8_2577,displ_4MN_1em8_2577,'*-k','LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (mm)')
legend('100kN: scale factor=40','4MN')
axis([0 1 -50 50])
text(0.8,30,'k^{hat} = 10^{-8} m^2/Pa.s','FontSize',16)
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('impulse_100kN_1em5_2577.txt');
time_100kN_1em5_2577=sim(:,1);
displ_100kN_1em5_2577=sim(:,4)*1000*40;

sim=load('impulse_4MN_1em5_2577.txt');
time_4MN_1em5_2577=sim(:,1);
displ_4MN_1em5_2577=sim(:,4)*1000;

% node displ time
figure(6)
plot(time_100kN_1em5_2577,displ_100kN_1em5_2577,'.-b',time_4MN_1em5_2577,displ_4MN_1em5_2577,'*-r','LineWidth',1)
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(mm)')
legend('100kN:Scale Factor=40','4MN')
axis([0 1 -60 60])
text(0.8,50,'khat=1e-5 m2/(Pa.s)','FontSize',14)
set(gca,'FontName','Helvetica','FontSize',14)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('impulse_4MN_1em5_k_2577.txt');
time_4MN_1em5_k_2577=sim(:,1);
displ_4MN_1em5_k_2577=sim(:,4)*1000;

% node displ time
figure(7)
plot(time_4MN_1em5_2577,displ_4MN_1em5_2577,'.-b',time_4MN_1em5_k_2577,displ_4MN_1em5_k_2577,'*-r','LineWidth',1)
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(mm)')
legend('khat=1e-5 m2/(Pa.s)','k(n^f)')
axis([0 1 -60 60])
text(0.8,40,'NODE A','FontSize',14)
set(gca,'FontName','Helvetica','FontSize',14)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('impulse_4MN_1em5_2457.txt');
time_4MN_1em5_2457=sim(:,1);
displ_4MN_1em5_2457=sim(:,4)*1000;

sim=load('impulse_4MN_1em5_k_2457.txt');
time_4MN_1em5_k_2457=sim(:,1);
displ_4MN_1em5_k_2457=sim(:,4)*1000;

% node displ time
figure(8)
plot(time_4MN_1em5_2457,displ_4MN_1em5_2457,'.-b',time_4MN_1em5_k_2457,displ_4MN_1em5_k_2457,'*-r','LineWidth',1)
xlabel('TIME(seconds)')
ylabel('VERTICAL DISPLACEMENT(mm)')
legend('khat=1e-5 m2/(Pa.s)','k(n^f)')
axis([0 1 -25 10])
text(0.8,-15,'NODE B','FontSize',14)
set(gca,'FontName','Helvetica','FontSize',14)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

