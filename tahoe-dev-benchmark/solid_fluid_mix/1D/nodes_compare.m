
%close all
clear all

format short e

sim=load('compress_40kpa_nd_189.txt');
time_40=sim(:,1);
displ_40=sim(:,4);

%small strain analytical
time_sol40=[0 1];
displ_sol40=[-0.0093 -0.0093];

% node displ time
figure(1)
plot(time_40,displ_40,'.-k',time_sol40,displ_sol40,'--k','LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
legend('finite strain simulation','analytical steady-state solution (small strain)')
text(0.6,-0.0088,'f_0=40kPa')
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('compress_2MPa_nd_189.txt');
time_2=sim(:,1);
displ_2=sim(:,4);

sim=load('compress_2MPa_nd_k_189.txt');
time_2_k=sim(:,1);
displ_2_k=sim(:,4);

%small strain analytical
time_sol2=[0 1];
displ_sol2=[-0.465 -0.465];

sim=load('compress_4MPa_nd_189.txt');
time_4=sim(:,1);
displ_4=sim(:,4);

sim=load('compress_4MPa_nd_k_189.txt');
time_4_k=sim(:,1);
displ_4_k=sim(:,4);

%small strain analytical
time_sol4=[0 1];
displ_sol4=[-0.93 -0.93];

sim=load('compress_8MPa_nd_189.txt');
time_8=sim(:,1);
displ_8=sim(:,4);

sim=load('compress_8MPa_nd_k_189.txt');
time_8_k=sim(:,1);
displ_8_k=sim(:,4);

%small strain analytical
time_sol8=[0 1];
displ_sol8=[-1.86 -1.86];

% node displ time
figure(2)
plot(time_2_k,displ_2_k,'o-k',time_2,displ_2,'.-k',time_sol2,displ_sol2,'--k',...
     time_4,displ_4,'.-k',time_4_k,displ_4_k,'o-k',time_sol4,displ_sol4,'--k',...
     time_8,displ_8,'.-k',time_8_k,displ_8_k,'o-k',time_sol8,displ_sol8,'--k',...
    'LineWidth',1)
xlabel('TIME (sec)')
ylabel('VERTICAL DISPLACEMENT (m)')
legend('k(n^f)','k^{hat} = const')
text(0.6,-0.38,'f_0=2MPa')
text(0.6,-0.75,'f_0=4MPa')
text(0.6,-1.45,'f_0=8MPa')
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})
