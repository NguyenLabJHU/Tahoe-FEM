% $Id: createValStruct.m,v 1.2 2001-08-14 21:41:54 strohban Exp $

function nDat = createValstruct(data,valDescr)
% nDat = createValStruct(nodaldata,valDescr)
% creates the nodal structure with all values
% nDat(timestep,nodenumber).D_X
% nDat(timestep,nodenumber).s11 etc.
  
  for step=1:size(data,1)
    fprintf(1,'Step %d (out of %d)\n',step,size(data,1));
    for outgroup=1:size(data,2)
      fprintf(1,'  Group Number %d (out of %d)\n',outgroup,size(data,2));
      d=data{step,outgroup};
      des=str2array(valDescr{outgroup});
      for j=1:length(des)
	if strcmp(des{j},'node')
	  % do nothing
	elseif strcmp(des{j},'D_X')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).D_X=d(i,j);
	  end 
	elseif strcmp(des{j},'D_Y')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).D_Y=d(i,j);
	  end
	elseif strcmp(des{j},'s11')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).s11=d(i,j);
	  end
	elseif strcmp(des{j},'s22')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).s22=d(i,j);
	  end
	elseif strcmp(des{j},'s12')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).s12=d(i,j);
	  end
	elseif strcmp(des{j},'d_t')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).d_t=d(i,j);
	  end
	elseif strcmp(des{j},'d_n')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).d_n=d(i,j);
	  end
	elseif strcmp(des{j},'T_t')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).T_t=d(i,j);
	  end
	elseif strcmp(des{j},'T_n')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).T_n=d(i,j);
	  end
	elseif strcmp(des{j},'lambda')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).lambda=d(i,j);
	  end
	elseif strcmp(des{j},'alpha')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).alpha=d(i,j);
	  end
	elseif strcmp(des{j},'norm_beta')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).norm_beta=d(i,j);
	  end
	elseif strcmp(des{j},'VM_Kirch')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).VM_Kirch=d(i,j);
	  end
	elseif strcmp(des{j},'press')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).press=d(i,j);
	  end
	elseif strcmp(des{j},'jump')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).jump=d(i,j);
	  end
	elseif strcmp(des{j},'Tmag')
	  for i=1:size(d,1)
	    nDat(step,d(i,1)).Tmag=d(i,j);
	  end
	else
	  fprintf(1,'Ignored: %s\n',des{j});
	end
      end
    end
  end
  
  





