% $Id: plotMesh.m,v 1.2 2001-08-14 21:41:54 strohban Exp $

function plotMesh (xy, elementsets, setnr, intColor,extColor,nodeNumbers,dotted);
% plotMesh (xy, elementsets, setnr, intColor,extColor,nodeNumbers,dotted)
% plots the meshes
% if setnr = [a,b,c] it plots sets a b and c.
% if setnr = [a:1:b] it plots sets from a to b.
% if setnr = 0 it plots all sets.
% inColor, extColor can be a color or "cycle" which will use
% a different color for each element set
% intColor= 'trans' makes it transparent
% if nodeNumbers = 1 it also displays the node numbers
% if dotted = 1 it plots dotted outlines

grp=elementsets;

numGroups=length(grp);
c = 'brgcmky';
numColors = length(c);

currColIndex=0;

range = setnr;

if setnr == 0
 range=[1:1:numGroups];
end

for g=range

fprintf(1,'Plotting Element Set %d',g);
cInt=intColor;
cExt=extColor;

if strcmp('cycle',intColor) | strcmp('cycle',extColor)
  currColIndex = mod(currColIndex+1,numColors);
  if strcmp('cycle',extColor)
    cExt=c(currColIndex+1);
  end
  if strcmp('cycle',intColor)
    cInt=c(currColIndex+1);
  end
end


for i = 1:size(grp{g},1)
   x = xy(grp{g}(i,:),1);
   y = xy(grp{g}(i,:),2);
   if strcmp(intColor,'trans')
   h = patch(x,y,'w');
   set(h,'FaceColor','none');   
   else
   h = patch(x,y,cInt);
   set(h,'FaceColor','none'); 
   end
   set(h,'EdgeColor',cExt);
   set(h,'LineWidth',1);
   if dotted
   set(h,'LineStyle','--');
   end

if nodeNumbers   
   for j=1:size(x,1)
     numX=x(j);
     numY=y(j);
     nodeNum=grp{g}(i,j);
     text(numX,numY,num2str(nodeNum),'FontSize',12);
end
   end
end
fprintf(1,' done \n');

end













