% $Id: readGeom.m,v 1.2 2001-08-14 21:41:54 strohban Exp $

function [xy,elementsets,nodesets] = readGeom(fname)
% [xy,elementsets,nodesets] = readGeom(fname) 
% returns nodal data and element data from new
% tahoe geometry input file. (usually ends in .geom)
  
  fid = fopen(fname,'rt');
  
  forwardUntil(fid,'*dimensions');
  numNodes=readInt(fid);
  numSpatialDims=readInt(fid);
  numElSets=readInt(fid);
  for i=1:numElSets+1
    line=readLine(fid); % dummy
  end
  numNodeSets=readInt(fid);
  
  fprintf(1,'Reading Nodesets (%d total)\n',numNodeSets);
  forwardUntil(fid,'*nodesets'); 
  for i=1:numNodeSets
    fprintf(1,'   Node Set: %d (out of %d)  ',i,numNodeSets);
    nodesets{i}=fetchNodeSet(fid);
    fprintf(1,'done\n');
  end
  
  fprintf(1,'Reading Elementsets (%d total)\n',numElSets);
  forwardUntil(fid,'*elements'); 
  for i=1:numElSets
    fprintf(1,'   Element Set: %d (out of %d)  ',i,numElSets);
    elementsets{i}=fetchElSet(fid);
    fprintf(1,'done\n');
  end
  
  fprintf(1,'Reading Nodes\n');
  xy=fetchNodes(fid);
  fprintf(1,'done\n');
  
  
function xy=fetchNodes(fid);
  forwardUntil(fid,'*nodes');
  numNodes=readInt(fid);
  numSpatDimsNodes=readInt(fid);
  xy=zeros(numNodes,numSpatDimsNodes);
  for i=1:numNodes
    if (numSpatDimsNodes == 2)
      line=readLine(fid);
      [data,count]=sscanf(line,'%d %f %f',3);
      if count ~= 3
	fprintf(1,'Problem reading Nodes');
	break;
      end
      xy(data(1),:) = data(2:3)';
    end
    if (numSpatDimsNodes == 3)
      line=readLine(fid);
      [data,count]=sscanf(line,'%d %f %f %f',4);
      if count ~= 4
	fprintf(1,'Problem reading Nodes');
	break;;
      end
      xy(data(1),:) = data(2:4)';
    end
  end
  
  
function group=fetchElSet(fid)
  forwardUntil(fid,'*set');
  numElsInSet=readInt(fid);
  numElNodes=readInt(fid);
  for i=1:numElsInSet
    line=readLine(fid);
    [data,count]=sscanf(line,'%d',numElNodes+1);
    if count ~= numElNodes+1
      fprintf(1,'Problem reading Element Set.');
      break;
    end
    group(i,:) = data(2:numElNodes+1)';
  end
  
function nodeSet=fetchNodeSet(fid)
  forwardUntil(fid,'*set');
  numNodes=readInt(fid);
  c=0; d=[];
  while c < numNodes
    line=readLine(fid);
    [data,count]=sscanf(line,'%d',inf);
    c=c+count;
    d=[d; data];
  end
  if (c ~= numNodes)
    fprintf(1,'Problem reading node set\n');
    return;
  end
  nodeSet = d';
  
function forwardUntil(fid,lookFor)
  while 1
    line = readLine(fid);
    if strcmp(lookFor,line)
      break;
    end
  end
  
function myInt = readInt(fid)
  line = readLine(fid);
  [myInt,count]=sscanf(line,'%d',inf);
  if count <= 0
    fprintf(1,'Problem reading Integer with line %s\n',line);
    break;
  end
  
function ll = readLine(fid)
  while 1
    ll = fgetl(fid);
    
    if ll < 0
      fprintf(1,'Error in readLine');
      return;
    end
    
    % remove all comments from line
    loc=findstr(ll,'#');
    if ~isempty(loc)
      ll=ll(1:loc-1);
    end
    
    % remove trailing blanks
    ll = deblank(ll);
    if ~isempty(ll)
      return;
    end
  end

  
  
  
  












