%close all
clc
clear all

format short e
% 
 Kdd=zeros(81,81);
 Kdphi=zeros(81,72);
 Kphid=zeros(72,81);
 Kphiphi=zeros(72,72);

 fK=zeros(153,153);
 
sim=load('fs_micromorph3D.info');

count=1;
for i=1:81
   for j=1:81
        
        Kdd(i,j)=sim(count);
        count=count+1;
    end
end
%%%%%%%%%%%%%
for i=1:81
    for j=1:72
        
        Kdphi(i,j)=sim(count);
        count=count+1;
    end
end
%%%%%%%%%%%%%%%  
for i=1:72
     for j=1:81   
        Kphid(i,j)=sim(count);
        count=count+1;
    end
end
% %%%%%%%%%%%%%%%%
for i=1:72
    for j=1:72
        Kphiphi(i,j)=sim(count);
        count=count+1;
    end
end
% 
% 
 %fK=[Kdd Kdphi; Kphid  Kphiphi];
%  fK=Kdd;
%  V=eig(fK);
%  RV=real(V);
%  
%  mv=max(V);
%  nV=V.*1/mv;
%  RnV=real(nV);
 
 BC1=[2 5 14 17 26 38 50 53 68 1 10 13 22 34 46 49 58 73 3 6 9 12 27 30 33 36 63];
BC=[1 2 3 5 6 9 10 12 13 14 17 22 26 27 30 33 34 36 38 46 49 50 53 58 63 68 73];
%BC=[3 6 9 12 27 30 33 36 63];


fK=[Kdd Kdphi; Kphid  Kphiphi];
length(BC);
n=1;
KK=fK;
while n<=length(BC)
   del=BC(n);
   del=del-(n-1);
   KK(:,del)=[];
   KK(del,:)=[];
   n=n+1;
        
end
% 
 %V=eig(KK);
 [V,D]=eig(KK);
 reV=real(V);
  rank(KK);
 
 for i=55:126
     
     for j=55:126
         
         phishape(i-54,j-54)=reV(i,j);
     end
 end



