
%close all
clear all

format short e

sim=load('nddataout_6_drained.txt');
time_drained=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3_drained=sim(:,5);

sim=load('nddataout_6_consol.txt');
time_consol=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3_consol=sim(:,5);

sim=load('nddataout_6_consol_elas.txt');
time_consol_elas=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3_consol_elas=sim(:,5);

sim=load('nddataout_6_dyn_step.txt');
time_dyn_step=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3_dyn_step=sim(:,5);

sim=load('nddataout_6_dyn_pi.txt');
time_dyn_pi=sim(:,1);
porepress=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3_dyn_pi=sim(:,5);

% node displ time
figure(1)
%plot(time_drained,d3_drained,'-k',time_consol,d3_consol,'--k',time_dyn,d3_dyn,':k','LineWidth',2)
%plot(time_consol_elas,d3_consol_elas,'-k',time_consol,d3_consol,'--k',time_dyn,d3_dyn,':k','LineWidth',2)
plot(time_consol_elas,d3_consol_elas,'-k',time_consol,d3_consol,'--k',time_dyn_step,d3_dyn_step,':k',time_dyn_pi,d3_dyn_pi,'-.k','LineWidth',2)
xlabel('time')
ylabel('d3')
legend('elastic consolidation','plastic consolidation','plastic dynamic step','plastic dynamic \pi')
set(gca,'FontName','Helvetica','FontSize',16)
set(gca,'YTickLabel',{'-10','-8','-6','-4','-2','0','2'})

sim=load('nddataout_2_drained.txt');
time_drained=sim(:,1);
porepress_drained=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_2_consol.txt');
time_consol=sim(:,1);
porepress_consol=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_2_consol_elas.txt');
time_consol_elas=sim(:,1);
porepress_consol_elas=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_2_dyn_step.txt');
time_dyn_step=sim(:,1);
porepress_dyn_step=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

sim=load('nddataout_2_dyn_pi.txt');
time_dyn_pi=sim(:,1);
porepress_dyn_pi=(sim(:,2))/1000;
d1=sim(:,3);
d2=sim(:,4);
d3=sim(:,5);

% node pore pressure time
figure(2)
%plot(time_drained,porepress_drained,'-k',time_consol,porepress_consol,'--k',time_dyn,porepress_dyn,':k','LineWidth',2)
%plot(time_consol_elas,porepress_consol_elas,'-k',time_consol,porepress_consol,'--k',time_dyn,porepress_dyn,':k','LineWidth',2)
plot(time_consol_elas,porepress_consol_elas,'-k',time_consol,porepress_consol,'--k',time_dyn_step,porepress_dyn_step,':k',time_dyn_pi,porepress_dyn_pi,'-.k','LineWidth',2)
xlabel('time')
ylabel('pore pressure')
%legend('drained','consolidating','dynamic')
%legend('elastic consolidation','plastic consolidation','plastic dynamic')
legend('elastic consolidation','plastic consolidation','plastic dynamic step','plastic dynamic pi')
set(gca,'FontName','Helvetica','FontSize',16)


% schedule function
figure(3)
time1=[0 1 2 2.5 3];
func1=[0 1 1 1   1];
time2=[0 1 2 2.5 3];
func2=[0 1 1 2   2];
time_dyn=0:0.01:3;
func_dyn=0.5 - 0.5*cos(pi*time_dyn)+time_dyn/3;
%func_dyn=0.5 - 0.5*cos(2*pi*time_dyn)+time_dyn/3;
%func_dyn=0.5 - 0.5*cos(3*pi*time_dyn)+time_dyn/3;
plot(time1,func1,':k',time2,func2,'-k',time_dyn,func_dyn,'--k','LineWidth',2)
xlabel('time')
ylabel('schedule function')
legend('step1','step2','\omega=\pi')
axis=([0 3 0 2.5])
set(gca,'FontName','Helvetica','FontSize',16)

