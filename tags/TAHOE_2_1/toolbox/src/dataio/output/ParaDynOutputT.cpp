#include "ParaDynOutputT.h"

#include "ParaDynT.h"
#include "OutputSetT.h"
#include "AutoArrayT.h"
#include <fstream.h>

#include "StringT.h"
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

ParaDynOutputT::ParaDynOutputT (ostream& out, 
				const ArrayT<StringT>& out_strings) :
  OutputBaseT (out, out_strings)
{
}

void ParaDynOutputT::WriteGeometry (void)
{
  ParaDynT par (fout);
  
  // write geometry file
  ofstream geo;
  OpenGeometryFile (par, geo);
  if (geo)
    {
      geo << fCoordinates->MajorDim() << "\n";
         
     for (int i=0; i < fElementSets.Length(); i++)
	 WritePart (geo,par,i);
    }
}

// *************** PRIVATE ********************

void ParaDynOutputT::OpenGeometryFile (ParaDynT& par, ofstream& geo) const
{
  StringT label = "atoms";
  StringT geofile;

  geofile = CreateFileName (label);
  geo.open (geofile);
  
  // header
  int h = 0;
  ArrayT<StringT> header;

  header.Allocate (1);
  header[h] = " ITEM: NUMBER OF ATOMS";

  par.WriteHeader (geo, header);
}

StringT ParaDynOutputT::CreateFileName (const StringT& Label) const
{
  StringT var (fOutroot);
  
  /* tack on extension */
  var.Append (".");
  var.Append (Label);
  return var;
}

void ParaDynOutputT::WriteBounds (ostream& geo, const ParaDynT& par) const
{
  dArray2DT local (fBounds->MajorDim(), fBounds->MinorDim());
  for (int i=0; i < local.MajorDim(); i++)
    for (int j=0; j < local.MinorDim(); j++)
      local(i,j) = (*fBounds)(i,j);

  par.WriteBoundHeader(geo);
  par.WriteBounds(geo,local);
}


void ParaDynOutputT::WritePart (ostream& geo, ParaDynT& par,int index) const
{
  WriteBounds(geo, par);
  par.WriteTime(geo);

  const ArrayT<StringT>& blockIDs = fElementSets[index]->BlockID();
  for (int b=0; b < fElementSets[index]->NumBlocks(); b++)
    {
      const iArrayT& nodes_used = fElementSets[index]->BlockNodesUsed(blockIDs[b]);
      WriteCoordinates (geo, par, nodes_used);
    }
}

void ParaDynOutputT::WriteCoordinates (ostream& geo, ParaDynT& par, 
				       const iArrayT& nodes) const
{
  dArray2DT local (nodes.Length(), fCoordinates->MinorDim());
  for (int i=0; i < nodes.Length(); i++)
    local.SetRow (i, (*fCoordinates)(nodes[i]));

  
  iArrayT tmp (fTypes->Length());
  for (int i=0; i < fTypes->Length(); i++)
    tmp[i] = (*fTypes)[i];

  iArrayT tmp2 (fParts->Length());
  for (int i=0; i < fParts->Length(); i++)
    tmp2[i] = (*fParts)[i];

  par.WriteCoordinateHeader (geo);
  par.WriteCoordinates (geo, local, tmp,tmp2);
}


