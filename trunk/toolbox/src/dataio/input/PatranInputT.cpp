/* $Id: PatranInputT.cpp,v 1.14 2003-11-10 22:14:22 cjkimme Exp $ */
/* created: sawimme July 2001 */

#include "PatranInputT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"


using namespace Tahoe;

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
  if (!fPatran.NamedComponents (names)) throw ExceptionT::kDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      if (!fPatran.ReadElementBlockDims (names[i], numelems, numelemnodes)) 
	throw ExceptionT::kDatabaseFail;
      if (numelems > 0)
	  groupnames[count++] = names[i];
    }
}

void PatranInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
  int count = 0;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  if (!fPatran.NamedComponents (names)) throw ExceptionT::kDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      int num = 0;
      if (!fPatran.NumNodesInSet (names[i], num)) throw ExceptionT::kDatabaseFail;
      if (num > 0)
	nodenames[count++] = names[i];
    }
}

int  PatranInputT::NumElementGroups (void) const
{
  int count = 0, numelems, numelemnodes;
  int numcomps = fPatran.NumNamedComponents ();
  ArrayT<StringT> names (numcomps);
  if (!fPatran.NamedComponents (names)) throw ExceptionT::kDatabaseFail;
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
  if (!fPatran.NamedComponents (names)) throw ExceptionT::kDatabaseFail;
  for (int i=0; i < numcomps; i++)
    {
      int num = 0;
      if (!fPatran.NumNodesInSet (names[i], num)) throw ExceptionT::kDatabaseFail;
      if (num > 0)
	count++;
    }
  return count;
}

void PatranInputT::ReadNodeID(iArrayT& node_id)
{
  if (!fPatran.ReadGlobalNodeMap(node_id)) throw ExceptionT::kDatabaseFail;
}

void PatranInputT::ReadCoordinates (dArray2DT& coords)
{
  if (!fPatran.ReadCoordinates (coords, coords.MinorDim()))
    throw ExceptionT::kDatabaseFail;
}

void PatranInputT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
  ReadCoordinates(coords);
  ReadNodeID(node_id);
}

int PatranInputT::NumElements (const StringT& name)
{
  int num, numnodes;
  if (!fPatran.ReadElementBlockDims (name, num, numnodes))
    throw ExceptionT::kDatabaseFail;
  return num;
}

int PatranInputT::NumElementNodes (const StringT& name)
{
  int num, numnodes;
  if (!fPatran.ReadElementBlockDims (name, num, numnodes))
    throw ExceptionT::kDatabaseFail;
  return numnodes;  
}

void PatranInputT::ReadAllElementMap (iArrayT& elemmap)
{
  cout << "\n\n PatranInputT::Not programmed to read all element map\n\n";
  elemmap = -1;
}

void PatranInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  PatranT::NamedTypes namedtype;
  if (!fPatran.ReadElementSet (name, namedtype, elemmap))
    throw ExceptionT::kDatabaseFail;
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
      if (index < 0 || index >= map.Length()) throw ExceptionT::kOutOfRange;
      set[n] = index;
    }  
}

void PatranInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  PatranT::NamedTypes namedtype;
  if (!fPatran.ReadConnectivity (name, namedtype, connects))
    throw ExceptionT::kDatabaseFail;

  /* convert from discontinuous to continuous numbering */
  iArrayT map (NumNodes());
  ReadNodeID(map);

  int *pc = connects.Pointer();
  for (int i=0; i < connects.Length(); i++, pc++)
    {
      int kdex;
      map.HasValue (*pc, kdex);
      if (kdex < 0 || kdex >= map.Length()) throw ExceptionT::kOutOfRange;
      *pc = kdex;
    }
}

void PatranInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& code)
{
  iArrayT elems;
  PatranT::NamedTypes namedtype;
  if (!fPatran.ReadElementSet (name, namedtype, elems))
    throw ExceptionT::kDatabaseFail;

  SetCode (namedtype, code);
}


int PatranInputT::NumNodesInSet (const StringT& name)
{
  int num;
  if (!fPatran.NumNodesInSet (name, num)) throw ExceptionT::kDatabaseFail;
  return num;
}

void PatranInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
  if (!fPatran.ReadNodeSet (name, nodes)) throw ExceptionT::kDatabaseFail;

  // offset and map to start numbering at zero
  // account for discontinuous numbering
  iArrayT map;
  ReadNodeID(map);
  for (int n=0; n < nodes.Length(); n++)
    {
      int index;
      map.HasValue (nodes[n], index);
      if (index < 0 || index >= map.Length()) throw ExceptionT::kOutOfRange;
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
  throw ExceptionT::kDatabaseFail;
}

void PatranInputT::ReadSideSetGlobal (const StringT& name, iArray2DT& sides) const
{
#pragma unused(name)
#pragma unused(sides)

  cout << "\n\n PatranInputT::Not programmed to read side sets\n\n";
  throw ExceptionT::kDatabaseFail;
}

/**************** PRIVATE *******************/

void PatranInputT::SetCode (PatranT::NamedTypes namedtype, GeometryT::CodeT& code) const
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
