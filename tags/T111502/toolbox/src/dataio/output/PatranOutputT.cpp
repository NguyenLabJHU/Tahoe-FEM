/* created: sawimme April 2002 */

#include "PatranOutputT.h"
#include "ofstreamT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "dArrayT.h"

using namespace Tahoe;

PatranOutputT::PatranOutputT (ostream& out, const ArrayT<StringT>& out_strings, bool binary) :
  OutputBaseT (out, out_strings),
  fBinary (binary),
  fPatran (fout)
{
  /* binary format not implemented yet */
  fBinary = false;
}

void PatranOutputT::WriteGeometry (void)
{
  StringT file;
  FileName (0, file, ".out");
  ofstreamT out (file);

  int num_elems = 0;
  int num_sets = fElementSets.Length();
  for (int e=0; e < num_sets; e++)
    num_elems += fElementSets[e]->NumElements();

  fPatran.WriteHeader (out, fCoordinates->MajorDim(), num_elems, fTitle);
  int firstID = 1;
  fPatran.WriteCoordinates (out, *fCoordinates, firstID);

  // write elements
  firstID = 1;
  iArrayT nodes_used (fCoordinates->MajorDim());
  nodes_used.SetValueToPosition();
  for (int e=0; e < num_sets; e++)
    if (fElementSets[e]->NumNodes() > 0)
      WriteConnectivity (out, firstID, e, nodes_used);

  // write named components
  CreateElementBlockIDs ();
  firstID = 1;
  for (int e=0; e < num_sets; e++)
    if (fElementSets[e]->NumNodes() > 0)
      WriteNamedComponents (out, firstID, e);

  fPatran.WriteClosure(out);
}

void PatranOutputT::WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
  OutputBaseT::WriteOutput (time, ID, n_values, e_values);
  StringT patfile;

  if (fElementSets[ID]->NumNodes() > 0)
    {
      CreateElementBlockIDs ();

      // open geometry file
      FileName (ID, patfile, ".out");
      ofstreamT patout (patfile);
      int numnodes = fElementSets[ID]->NumNodes();
      int numelems = fElementSets[ID]->NumElements();
      fPatran.WriteHeader (patout, numnodes, numelems, fTitle);

      // write nodes
      int firstID = 1;
      iArrayT nodes_used;
      nodes_used.Alias (fElementSets[ID]->NodesUsed());
      dArray2DT coords (nodes_used.Length(), fCoordinates->MinorDim());
      coords.RowCollect (nodes_used, *fCoordinates);
      fPatran.WriteCoordinates (patout, coords, firstID);

      // write elements
      WriteConnectivity (patout, firstID, ID, nodes_used);

      // write named components
      firstID = 1;
      WriteNamedComponents (patout, firstID, ID);

      fPatran.WriteClosure(patout);
    }
}

/*************************************************************************
 * Private
 *************************************************************************/

/* generate database file name for the given ID */
void PatranOutputT::FileName (int ID, StringT& filename, const char* ext) const
{
  /* root */
  filename = fOutroot;
  
  /* tack on sequence number */
  if (fSequence > 0) filename.Append (".sq", fSequence + 1);

  /* I/O ID */
  filename.Append (".io", ID);

  /* changing geometry */
  if (fElementSets[ID]->Changing())
    filename.Append (".ps", fElementSets[ID]->PrintStep() + 1);

  /* extension */
  filename.Append (ext);
}

void PatranOutputT::WriteConnectivity (ostream& patout, int& firstID, int ID, iArrayT& nodes_used) const
{
  int num_blocks = fElementSets[ID]->NumBlocks();
  for (int i=0; i < num_blocks; i++)
    {
      StringT blockid = fElementSets[ID]->BlockID(i);
      ArrayT<PatranT::ElementTypes> types (fElementSets[ID]->NumElements());
      types = GetPatranElementType (fElementSets[ID]->Geometry());
      
      // write connectivity
      const iArray2DT* connects = fElementSets[ID]->Connectivities(blockid);
      iArray2DT localconn (connects->MajorDim(), connects->MinorDim());
      LocalConnectivity (nodes_used, *connects, localconn);
      localconn++;
      fPatran.WriteElements (patout, localconn, types, firstID);
      
      // increment element ID
      firstID += connects->MajorDim();
    }
}

void PatranOutputT::WriteNamedComponents (ostream& patout, int& firstID, int ID) const
{
  int num_blocks = fElementSets[ID]->NumBlocks();
  for (int i=0; i < num_blocks; i++)
    {
      int types = GetPatranElementType (fElementSets[ID]->Geometry());
      iArrayT eids (fElementSets[ID]->NumElements());
      eids.SetValueToPosition();
      eids += firstID;
      iArray2DT comps (fElementSets[ID]->NumElements(), 2);
      comps.SetColumn (0, types);
      comps.SetColumn (1, eids);
      
      int bid = fElementBlockIDs[ID][i];
      fPatran.WriteNamedComponent (patout, fElementSets[ID]->BlockID(i), bid, comps);
      
      // increment element ID
      firstID += fElementSets[ID]->NumElements();
    }  
}


PatranT::NamedTypes PatranOutputT::GetPatranNamedType (GeometryT::CodeT geom) const
{
  switch (geom)
    {
    case GeometryT::kLine:          return PatranT::kNCLine;
    case GeometryT::kQuadrilateral: return PatranT::kNCQuad;
    case GeometryT::kTriangle:      return PatranT::kNCTriangle;
    case GeometryT::kHexahedron:    return PatranT::kNCHex;
    case GeometryT::kTetrahedron:   return PatranT::kNCTet;
    case GeometryT::kPentahedron:   return PatranT::kNCWedge;
    }
  return PatranT::kNoNamedType;
}

PatranT::ElementTypes PatranOutputT::GetPatranElementType (GeometryT::CodeT geom) const
{
  switch (geom)
    {
    case GeometryT::kLine:          return PatranT::kLine;
    case GeometryT::kQuadrilateral: return PatranT::kQuadrilateral;
    case GeometryT::kTriangle:      return PatranT::kTriangle;
    case GeometryT::kHexahedron:    return PatranT::kHexahedron;
    case GeometryT::kTetrahedron:   return PatranT::kTetrahedron;
    case GeometryT::kPentahedron:   return PatranT::kPentahedron;
    }
  return PatranT::kNoElementType;
}
