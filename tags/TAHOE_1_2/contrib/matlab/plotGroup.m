% $Id: plotGroup.m,v 1.2 2001-08-14 21:41:54 strohban Exp $

function plotGroup (xy, elementsets,setnr, nDat, step,type,edgecolor)
% plotGroup (xy, elementsets,setnr, nDat, step,type, edgecolor)
% plots values for timestep step 
% the element set is selected by setnr 
% the values plotted depend on type
% i.e type = 's11' plots 11 stress etc.
% egdecolor = 0 does not print any edges
% edgecolor = 'k' prints black element outlines
  
  for i = 1:size(elementsets{setnr},1);
    x = xy(elementsets{setnr}(i,:),1)+[nDat(step,elementsets{setnr}(i,:)).D_X]';
    y = xy(elementsets{setnr}(i,:),2)+[nDat(step,elementsets{setnr}(i,:)).D_Y]';
    
    % get the appropriate values out of the structure
    if strcmp(type,'D_X')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).D_X]);      
    elseif strcmp(type,'D_Y')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).D_Y]);
    elseif strcmp(type,'s11')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).s11]);
    elseif strcmp(type,'s22')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).s22]);
    elseif strcmp(type,'s12')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).s12]);
    elseif strcmp(type,'d_t')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).d_t]);  
    elseif strcmp(type,'d_n')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).d_n]);
    elseif strcmp(type,'T_t')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).T_t]);
    elseif strcmp(type,'T_n')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).T_n]);
    elseif strcmp(type,'lambda')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).lambda]);
    elseif strcmp(type,'alpha')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).alpha]);
    elseif strcmp(type,'norm_beta')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).norm_beta]);
    elseif strcmp(type,'VM_Kirch')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).VM_Kirch]);
    elseif strcmp(type,'press')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).press]);
    elseif strcmp(type,'jump')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).jump]);
    elseif strcmp(type,'Tmag')
      val=transpose([nDat(step,elementsets{setnr}(i,:)).Tmag]);
    else
      fprintf(1,'unknown type %s\n',type);
      val='c';
    end
    
    % actually draw the patch
    h = patch(x,y,val);
    
    if edgecolor == 0
      set(h,'EdgeColor','none');
    else
      set(h,'EdgeColor',edgecolor);
      set(h,'LineWidth',1);
    end
  end
  
  

