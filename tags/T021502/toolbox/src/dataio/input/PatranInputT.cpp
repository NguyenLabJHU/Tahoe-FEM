/* $Id: PatranInputT.cpp,v 1.8 2002-01-27 18:38:12 paklein Exp $ */
/* created: sawimme July 2001 */

#include "PatranInputT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

PatranInputT::PatranInputT (ostream& out) :
  InputBaseT (out),
  fPatran (out)
{
}

bool PatranInputT::Open (const StringT& file)
{
	if (!fPatran.OpenRead (file)) {
		cout << "\n PatranInputT::Open: error opening file: " << file << endl;
		return false;
	} else return true;
}

void PatranInputT::Close (void)
{
}

void PatranInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
  int count = 0, numelems, numelemnodes;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  if (!fPatran.NamedComponents (names)) throw eDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadElementBlockDims (names[i], numelems, numelemnodes)) 
	throw eDatabaseFail;
      if (numelems > 0)
	  groupnames[count++] = names[i];
    }
}

void PatranInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
  int count = 0;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  iArrayT nodes;
  if (!fPatran.NamedComponents (names)) throw eDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadNodeSet (names[i], nodes)) throw eDatabaseFail;
      if (nodes.Length() > 0)
	nodenames[count++] = names[i];
    }
}

int  PatranInputT::NumElementGroups (void) const
{
  int count = 0, numelems, numelemnodes;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  if (!fPatran.NamedComponents (names)) throw eDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadElementBlockDims (names[i], numelems, numelemnodes)) 
	return false;
      if (numelems > 0)
	count++;
    }
  return count;
}

int  PatranInputT::NumNodeSets (void) const
{
  int count = 0;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  iArrayT nodes;
  if (!fPatran.NamedComponents (names)) throw eDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadNodeSet (names[i], nodes)) throw eDatabaseFail;
      if (nodes.Length() > 0)
	count++;
    }
  return count;
}

void PatranInputT::ReadNodeMap (iArrayT& nodemap)
{
  if (!fPatran.ReadGlobalNodeMap (nodemap)) throw eDatabaseFail;
}

void PatranInputT::ReadCoordinates (dArray2DT& coords)
{
  if (!fPatran.ReadCoordinates (coords, coords.MinorDim()))
    throw eDatabaseFail;
}

void PatranInputT::ReadCoordinates (dArray2DT& coords, iArrayT& nodemap)
{
  ReadCoordinates (coords);
  ReadNodeMap (nodemap);
}

int PatranInputT::NumElements (const StringT& name)
{
  int num, numnodes;
  if (!fPatran.ReadElementBlockDims (name, num, numnodes))
    throw eDatabaseFail;
  return num;
}

int PatranInputT::NumElementNodes (const StringT& name)
{
  int num, numnodes;
  if (!fPatran.ReadElementBlockDims (name, num, numnodes))
    throw eDatabaseFail;
  return numnodes;  
}

void PatranInputT::ReadAllElementMap (iArrayT& elemmap)
{
  cout << "\n\n PatranInputT::Not programmed to read all element map\n\n";
  elemmap = -1;
}

void PatranInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  int namedtype;
  if (!fPatran.ReadElementSet (name, namedtype, elemmap))
    throw eDatabaseFail;
}

void PatranInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
  ReadGlobalElementMap (name, set);

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map;
  ReadAllElementMap (map);
  for (int n=0; n < set.Length(); n++)
    {
      int index;
      map.HasValue (set[n], index);
      if (index < 0 || index >= map.Length()) throw eOutOfRange;
      set[n] = index;
    }  
}

void PatranInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  int namedtype;
  if (!fPatran.ReadConnectivity (name, namedtype, connects))
    throw eDatabaseFail;

  /* convert from discontinuous to continuous numbering */
  iArrayT map (NumNodes());
  ReadNodeMap (map);

  int *pc = connects.Pointer();
  for (int i=0; i < connects.Length(); i++, pc++)
    {
      int kdex;
      map.HasValue (*pc, kdex);
      if (kdex < 0 || kdex >= map.Length()) throw eOutOfRange;
      *pc = kdex;
    }
}

void PatranInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& code)
{
  iArrayT elems;
  int namedtype;
  if (!fPatran.ReadElementSet (name, namedtype, elems))
    throw eDatabaseFail;

  SetCode (namedtype, code);
}


int PatranInputT::NumNodesInSet (const StringT& name)
{
  int num;
  if (!fPatran.NumNodesInSet (name, num)) throw eDatabaseFail;
  return num;
}

void PatranInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
  if (!fPatran.ReadNodeSet (name, nodes)) throw eDatabaseFail;

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map;
  ReadNodeMap (map);
  for (int n=0; n < nodes.Length(); n++)
    {
      int index;
      map.HasValue (nodes[n], index);
      if (index < 0 || index >= map.Length()) throw eOutOfRange;
      nodes[n] = index;
    }
}

int PatranInputT::NumSidesInSet (const StringT& anme) const
{
#pragma unused(anme)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  return 0;
}

StringT PatranInputT::SideSetGroupName (const StringT& name) const
{
#pragma unused(name)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  StringT elname ("");
  return elname; 
}

void PatranInputT::ReadSideSetLocal (const StringT& name, iArray2DT& sides) const
{
#pragma unused(name)
#pragma unused(sides)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  throw eDatabaseFail;
}

void PatranInputT::ReadSideSetGlobal (const StringT& name, iArray2DT& sides) const
{
#pragma unused(name)
#pragma unused(sides)
  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  throw eDatabaseFail;
}

/**************** PRIVATE *******************/

void PatranInputT::SetCode (int namedtype, GeometryT::CodeT& code) const
{
  switch (namedtype)
    {
    case PatranT::kNCLine:
    case PatranT::kNCLine2:
    case PatranT::kNCLine3:
      {
	code = GeometryT::kLine;
	break;
      }
    case PatranT::kNCTriangle: 
    case PatranT::kNCTriangle2:
    case PatranT::kNCTriangle3:
      {
	code = GeometryT::kTriangle;
	break;
      }
    case PatranT::kNCQuad: 
    case PatranT::kNCQuad2:
    case PatranT::kNCQuad3:
       {
	 code = GeometryT::kQuadrilateral;
	 break;
       }
    case PatranT::kNCTet:
    case PatranT::kNCTet2:
    case PatranT::kNCTet3:
      {
	code = GeometryT::kTetrahedron;
	break;
      }
    case PatranT::kNCWedge: 
    case PatranT::kNCWedge2:
    case PatranT::kNCWedge3:
      {
	code = GeometryT::kPentahedron;
	break;
      }
    case PatranT::kNCHex:
    case PatranT::kNCHex2:
    case PatranT::kNCHex3:
      {
	code = GeometryT::kHexahedron;
	break;
      }
    default:
      code = GeometryT::kNone;
    }
}
