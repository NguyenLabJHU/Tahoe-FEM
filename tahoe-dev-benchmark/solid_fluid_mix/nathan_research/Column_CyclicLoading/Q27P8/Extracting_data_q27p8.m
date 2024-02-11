clc
clear

% Node on the Top %

Node_1 = load('result_1.txt');
S33_1 = Node_1(:,3);
p_f_1 = Node_1(:,4);
e33_1 = Node_1(:,5);

Node_13 = load('result_13.txt');
S33_13 = Node_13(:,3);
p_f_13 = Node_13(:,4);
e33_13 = Node_13(:,5);

Node_5 = load('result_5.txt');
S33_5 = Node_5(:,3);
p_f_5 = Node_5(:,4);
e33_5 = Node_5(:,5);

Node_260 = load('result_260.txt');
S33_260 = Node_260(:,3);
p_f_260 = Node_260(:,4);
e33_260 = Node_260(:,5);

Node_256 = load('result_256.txt');
S33_256 = Node_256(:,3);
p_f_256 = Node_256(:,4);
e33_256 = Node_256(:,5);

Node_12 = load('result_12.txt');
S33_12 = Node_12(:,3);
p_f_12 = Node_12(:,4);
e33_12 = Node_12(:,5);

Node_24 = load('result_24.txt');
S33_24 = Node_24(:,3);
p_f_24 = Node_24(:,4);
e33_24 = Node_24(:,5);

Node_20 = load('result_20.txt');
S33_20 = Node_20(:,3);
p_f_20 = Node_20(:,4);
e33_20 = Node_20(:,5);

Node_270 = load('result_270.txt');
S33_270 = Node_270(:,3);
p_f_270 = Node_270(:,4);
e33_270 = Node_270(:,5);

Node_267 = load('result_267.txt');
S33_267 = Node_267(:,3);
p_f_267 = Node_267(:,4);
e33_267 = Node_267(:,5);

Node_259 = load('result_259.txt');
S33_259 = Node_259(:,3);
p_f_259 = Node_259(:,4);
e33_259 = Node_259(:,5);

Node_263 = load('result_263.txt');
S33_263 = Node_263(:,3);
p_f_263 = Node_263(:,4);
e33_263 = Node_263(:,5);

Node_8 = load('result_8.txt');
S33_8 = Node_8(:,3);
p_f_8 = Node_8(:,4);
e33_8 = Node_8(:,5);

Node_16 = load('result_16.txt');
S33_16 = Node_16(:,3);
p_f_16 = Node_16(:,4);
e33_16 = Node_16(:,5);

Node_4 = load('result_4.txt');
S33_4 = Node_4(:,3);
p_f_4 = Node_4(:,4);
e33_4 = Node_4(:,5);

Node_364 = load('result_364.txt');
S33_364 = Node_364(:,3);
p_f_364 = Node_364(:,4);
e33_364 = Node_364(:,5);

Node_367 = load('result_367.txt');
S33_367 = Node_367(:,3);
p_f_367 = Node_367(:,4);
e33_367 = Node_367(:,5);

Node_165 = load('result_165.txt');
S33_165 = Node_165(:,3);
p_f_165 = Node_165(:,4);
e33_165 = Node_165(:,5);

Node_169 = load('result_169.txt');
S33_169 = Node_169(:,3);
p_f_169 = Node_169(:,4);
e33_169 = Node_169(:,5);

Node_160 = load('result_160.txt');
S33_160 = Node_160(:,3);
p_f_160 = Node_160(:,4);
e33_160 = Node_160(:,5);

Node_359 = load('result_359.txt');
S33_359 = Node_359(:,3);
p_f_359 = Node_359(:,4);
e33_359 = Node_359(:,5);

Node_361 = load('result_361.txt');
S33_361 = Node_361(:,3);
p_f_361 = Node_361(:,4);
e33_361 = Node_361(:,5);

Node_157 = load('result_157.txt');
S33_157 = Node_157(:,3);
p_f_157 = Node_157(:,4);
e33_157 = Node_157(:,5);

Node_162 = load('result_162.txt');
S33_162 = Node_162(:,3);
p_f_162 = Node_162(:,4);
e33_162 = Node_162(:,5);

Node_155 = load('result_155.txt');
S33_155 = Node_155(:,3);
p_f_155 = Node_155(:,4);
e33_155 = Node_155(:,5);

Node_82 = load('result_82.txt');
S33_82 = Node_82(:,3);
p_f_82 = Node_82(:,4);
e33_82 = Node_82(:,5);

Node_86 = load('result_86.txt');
S33_86 = Node_86(:,3);
p_f_86 = Node_86(:,4);
e33_86 = Node_86(:,5);

Node_64 = load('result_64.txt');
S33_64 = Node_64(:,3);
p_f_64 = Node_64(:,4);
e33_64 = Node_64(:,5);

Node_68 = load('result_68.txt');
S33_68 = Node_68(:,3);
p_f_68 = Node_68(:,4);
e33_68 = Node_68(:,5);

Node_46 = load('result_46.txt');
S33_46 = Node_46(:,3);
p_f_46 = Node_46(:,4);
e33_46 = Node_46(:,5);


Stress = load('Stress.txt');

S_33 = (S33_1+S33_13+S33_5+S33_260+S33_256+S33_257+S33_270+S33_20+S33_24+S33_12+...
    S33_259+S33_263+S33_8+S33_16+S33_4+S33_364+S33_367+S33_165+S33_169+S33_160+...
    S33_359+S33_361+S33_157+S33_162+S33_155)/25;
e_33 = (e33_1+e33_13+e33_5+e33_260+e33_256+e33_257+e33_270+e33_20+e33_24+e33_12+...
    e33_259+e33_263+e33_8+e33_16+e33_4+e33_364+e33_367+e33_165+e33_169+e33_160+...
    e33_359+e33_361+e33_157+e33_162+e33_155)/25;
p_f = (p_f_82+p_f_86+p_f_64+p_f_68+p_f_46)/5;



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


%%% The number 14 in the plot commands is related to the number of time
%%% step which covers the first cycle of loading function.
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



