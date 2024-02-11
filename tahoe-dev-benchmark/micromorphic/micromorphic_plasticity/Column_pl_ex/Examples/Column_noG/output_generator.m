
% number of Gauss points Nip
Nip=27;
% number of the output variables n
%n=43;
n=38; % with micro plasticity  parameters

%n=74; % 38 + 36 (F,Fe,X,Xe)=74
% number of y (YES) ny hits ny=Nip*n
ny=Nip*n;

% open the file 
output = fopen('out_Stan_Pl','w');
fprintf(output,'%% \n');
fprintf(output,'6 \n');
fprintf(output,'5 Column_ex_Standard_pl_NR.io0.exo \n');
fprintf(output,'0 \n');
fprintf(output,'el_Stan_Pl \n');
fprintf(output,'1 \n');

for i=1:ny
    fprintf(output,'y \n');
end

fprintf(output,'1 \n');
fprintf(output,'1 \n');

% number of elements for the results n_el
% if n_el =2 you should enter two el_num, if n_el=3, three el_num seperately 
n_el=2;
fprintf(output,'%d \n',n_el);

% element number el_num
% if chosen more than one element, each element should be entered manually
% first element number
el_num=1;
fprintf(output,'%d \n',el_num);
% second element number
% el_num=1; % 2,3,4,5,,,,,,721 etc.
% fprintf(output,'%d \n',el_num);
 el_num=8;
fprintf(output,'%d \n',el_num);   