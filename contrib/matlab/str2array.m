% $Id: str2array.m,v 1.2 2001-08-14 21:41:55 strohban Exp $

function A = str2array(str)
k=1;
c=1;
i=1;

aa=find(~isspace(str));
vold=aa(1);
while i<=length(aa)
     v=aa(i);
if v>vold+1
c=c+1;
k=1;
end
A{c}(k)=str(v);
k=k+1;
vold=v;
i=i+1;
end
