% $Id: readNodal.m,v 1.2 2001-08-14 21:41:54 strohban Exp $

function [t,nodaldata,valDescr] = readNodal(fname)
% [t,nodaldata,valDescr] = readNodal(fname)
% reads in nodal data from an .run file
% nodaldata contains the nodal data for each output set with indeces
% nodaldata{timestepnr,outputset}(row,column)
% t is an array with t(timestepnr) = time
% valDescr{groupnr} is a string explaining what each values are
% in the nodaldata array since different groups have different values
% and column names

fid = fopen(fname,'rt');

timestep = 0;

while 1
   line = fgetl(fid);
   if isempty(line) & feof(fid)
      break;
   end

   l = length(line);

   if l>= 5 & line(1:5)==' Time'
      time = sscanf(line,' Time. . . . . . . . . . . . . . . . . . . . . . = %f',1);
      if timestep == 0 | time > t(timestep)
         fprintf(1,'%s\n',line);
         timestep = timestep + 1;
         t(timestep) = time;
      end
   elseif l>=13 & line(2:13) == 'Group number'
      group = sscanf(line,...
                ' Group number. . . . . . . . . . . . . . . . . . = %d',1);
      fprintf(1,'%s\n',line);
      % fprintf(1,'   I think that means group %d\n',group);
   elseif l>= 23 & line(2:23) == 'Number of nodal points'
      numnodes = sscanf(line,...
                ' Number of nodal points. . . . . . . . . . . . . = %d',1);
      % fprintf(1,'%s\n',line);
      % fprintf(1,'   I think that means %d nodes\n',numnodes);
   elseif l>=17 & line(2:17) == 'Number of values'
      numvals = sscanf(line,...
                ' Number of values. . . . . . . . . . . . . . . . = %d',1);
      % fprintf(1,'%s\n',line);
      % fprintf(1,'   I think that means %d values\n',numvals);
   elseif l>=8 & line(5:8) == 'node'     
      valDescr{group}=line;
      [data,count] = fscanf(fid,'%f',[numvals+1 numnodes]);
       
      if count < (numvals+1)*numnodes
         fprintf(1,'   truncated file\n');
         nodaldata(timestep,:) = [];
         t(end) = [];
         return;
      end
      nodaldata{timestep,group} = data';
      
   end
end


fclose(fid);











