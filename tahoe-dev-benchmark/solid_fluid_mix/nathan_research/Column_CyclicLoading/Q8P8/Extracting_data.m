clc
clear

% Node on the Top %

Node_13 = load('result_13.txt');
S33_13 = Node_13(:,3);
p_f_13 = Node_13(:,4);
e33_13 = Node_13(:,5);

Node_17 = load('result_17.txt');
S33_17 = Node_17(:,3);
p_f_17 = Node_17(:,4);
e33_17 = Node_17(:,5);

Node_61 = load('result_61.txt');
S33_61 = Node_61(:,3);
p_f_61 = Node_61(:,4);
e33_61 = Node_61(:,5);

Node_64 = load('result_64.txt');
S33_64 = Node_64(:,3);
p_f_64 = Node_64(:,4);
e33_64 = Node_64(:,5);

Node_20 = load('result_20.txt');
S33_20 = Node_20(:,3);
p_f_20 = Node_20(:,4);
e33_20 = Node_20(:,5);

Node_16 = load('result_16.txt');
S33_16 = Node_16(:,3);
p_f_16 = Node_16(:,4);
e33_16 = Node_16(:,5);

Node_77 = load('result_77.txt');
S33_77 = Node_77(:,3);
p_f_77 = Node_77(:,4);
e33_77 = Node_77(:,5);

Node_49 = load('result_49.txt');
S33_49 = Node_49(:,3);
p_f_49 = Node_49(:,4);
e33_49 = Node_49(:,5);

Node_47 = load('result_47.txt');
S33_47 = Node_47(:,3);
p_f_47 = Node_47(:,4);
e33_47 = Node_47(:,5);

Node_25 = load('result_25.txt');
p_f_25 = Node_25(:,4);
time1 = Node_25(:,1);

Node_33 = load('result_33.txt');
p_f_33 = Node_33(:,4);

Node_29 = load('result_29.txt');
p_f_29 = Node_29(:,4);

Stress = load('Stress.txt');

S_33 = (S33_13+S33_17+S33_61+S33_16+S33_20+S33_64+S33_47+S33_49+S33_77)/9;
e_33 = (e33_13+e33_17+e33_61+e33_16+e33_20+e33_64+e33_47+e33_49+e33_77)/9;
p_f = (p_f_25+p_f_29+p_f_33)/3;



Exper = [25.1572	0.204778;
50.3145	0.511945;
75.4717	0.887372;
100.629	1.19454;
134.172	1.77474;
159.329	2.28669;
184.486	2.83276;
209.644	3.27645;
226.415	3.72014;
243.187	4.16382;
259.958	4.47099];



figure(1)
plot(time1(1:14,1),p_f(1:14,1),Exper(:,1),Exper(:,2)*1000,'g')
xlabel('Time (s)','fontsize',16)
ylabel('Pore water pressure (Pa)','fontsize',16)
grid on
legend('Tahoe mixture theory','Exprimental data')
set(gca,'FontName','Helvetica','FontSize',14)

figure(2)
plot(-e_33(1:14,1),-S_33(1:14,1),Stress(:,1)/100,Stress(:,2)*1000)
xlabel('Strain','fontsize',16)
ylabel('Stress (Pa)','fontsize',16)
grid on
legend('Tahoe mixture theory','Exprimental data')
set(gca,'FontName','Helvetica','FontSize',14)



