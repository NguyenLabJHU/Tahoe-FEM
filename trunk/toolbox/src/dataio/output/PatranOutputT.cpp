/* created: sawimme April 2002 */

#include "PatranOutputT.h"
#include "PatranT.h"
#include "ofstreamT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"

PatranOutputT::PatranOutputT (ostream& out, const ArrayT<StringT>& out_strings, bool binary) :
  OutputBaseT (out, out_strings),
  fBinary (binary)
{
  /* binary format not implemented yet */
  fBinary = false;
}

void PatranOutputT::WriteGeometry (void)
{
  PatranT pat (fout);
  StringT file;
  FileName (0, file, ".pao");
  ofstreamT out (file);

  int num_elems = 0;
  int num_sets = fElementSets.Length();
  for (int e=0; e < num_sets; e++)
    num_elems += fElementSets[e]->NumElements();

  pat.WriteHeader (out, fCoordinates->MajorDim(), num_elems, fTitle);
  int firstID = 1;
  pat.WriteCoordinates (out, *fCoordinates, firstID);

  // write elements
  for (int e=0; e < num_sets; e++)
    if (fElementSets[e]->NumNodes() > 0)
      {
	int num_blocks = fElementSets[e]->NumBlocks();
	for (int i=0; i < num_blocks; i++)
	  {
	    StringT blockid = fElementSets[e]->BlockID(i);
	    iArrayT types (fElementSets[e]->NumElements());
	    types = GetPatranElementType (fElementSets[e]->Geometry());
	    pat.WriteElements (out, *(fElementSets[e]->Connectivities (blockid)), types, firstID);
	    firstID += fElementSets[e]->NumElements();
	  }
      }

  // write named components
  firstID = 1;
  for (int e=0, id=0; e < num_sets; e++)
    if (fElementSets[e]->NumNodes() > 0)
      {
	int num_blocks = fElementSets[e]->NumBlocks();
	for (int i=0; i < num_blocks; i++, id++)
	  {
	    int types = GetPatranElementType (fElementSets[e]->Geometry());
	    iArrayT eids (fElementSets[e]->NumElements());
	    eids.SetValueToPosition();
	    eids += firstID;
	    iArray2DT comps (fElementSets[e]->NumElements(), 2);
	    comps.SetColumn (0, types);
	    comps.SetColumn (1, eids);
	    pat.WriteNamedComponent (out, fElementSets[e]->BlockID(i), id, comps);
	    firstID += fElementSets[e]->NumElements();
	  }
      }

  pat.WriteClosure(out);
}

void PatranOutputT::WriteOutput (double time, int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
  OutputBaseT::WriteOutput (time, ID, n_values, e_values);

  PatranT pat (fout);
  StringT patfile;

  if (fElementSets[ID]->NumNodes() > 0)
    {
      // write geometry file
      FileName (ID, patfile, "_geo.pao");
      ofstreamT patout (patfile);
      int numnodes = fElementSets[ID]->NumNodes();
      int numelems = fElementSets[ID]->NumElements();
      pat.WriteHeader (patout, numnodes, numelems, fTitle);

      // write nodes
      int firstID = 1;
      iArrayT nodes_used;
      nodes_used.Alias (fElementSets[ID]->NodesUsed());
      dArray2DT coords (nodes_used.Length(), fCoordinates->MinorDim());
      coords.RowCollect (nodes_used, *fCoordinates);
      pat.WriteCoordinates (patout, coords, firstID);

      // write elements
      int num_blocks = fElementSets[ID]->NumBlocks();
      for (int i=0; i < num_blocks; i++)
	{
	  StringT blockid = fElementSets[ID]->BlockID(i);
	  iArrayT types (fElementSets[ID]->NumElements());
	  types = GetPatranElementType (fElementSets[ID]->Geometry());

	  // write connectivity
	  const iArray2DT* connects = fElementSets[ID]->Connectivities(blockid);
	  iArray2DT localconn (connects->MajorDim(), connects->MinorDim());
	  LocalConnectivity (nodes_used, *connects, localconn);
	  localconn++;
	  pat.WriteElements (patout, localconn, types, firstID);

	  // write named components
	  iArrayT eids (fElementSets[ID]->NumElements());
	  eids.SetValueToPosition();
	  eids += firstID;
	  iArray2DT comps (fElementSets[ID]->NumElements(), 2);
	  comps.SetColumn (0, types);
	  comps.SetColumn (0, eids);
	  pat.WriteNamedComponent (patout, blockid, i+1, comps);

	  // increment element ID
	  firstID += connects->MajorDim();
	}
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

int PatranOutputT::GetPatranElementType (GeometryT::CodeT geom) const
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
  return GeometryT::kNone;
}
