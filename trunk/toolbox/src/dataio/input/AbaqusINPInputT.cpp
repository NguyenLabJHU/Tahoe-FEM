/* created: sawimme Oct 2006 */

#include "AbaqusINPInputT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"

/* Need to do yet
	1. NumGlobalElements
	2. ReadAllElementMap
	3. ReadGlobalElementSet
	4. NumSidesInSet
	5. SideSetGroupName
	6. ReadSideSetLocal
	7. ReadSideSetGlobal
*/


using namespace Tahoe;

AbaqusINPInputT::AbaqusINPInputT (ostream& out) :
  InputBaseT (out),
  fAbaqus ()
{
}

bool AbaqusINPInputT::Open (const StringT& file)
{
	if (!fAbaqus.OpenRead (file)) {
		cout << "\n AbaqusINPInputT::Open: error opening file: " << file << endl;
		return false;
	} else return true;
}

void AbaqusINPInputT::Close (void)
{
}

int AbaqusINPInputT::NumNodes (void) const
{
	return fAbaqus.NumNodes ();
}

void AbaqusINPInputT::ElementGroupNames (ArrayT<StringT>& groupnames) const
{
	fAbaqus.SetNames (groupnames, "ELSET", "ELEMENT");
}

void AbaqusINPInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{
	fAbaqus.SetNames (nodenames, "NSET", "NODE");
}

int AbaqusINPInputT::NumElementGroups (void) const
{
  return fAbaqus.NumElSets ();
}

int AbaqusINPInputT::NumNodeSets (void) const
{
  return fAbaqus.NumNodeSets ();
}

void AbaqusINPInputT::ReadNodeID(iArrayT& node_id)
{
	fAbaqus.NodeIDs(node_id);
}

void AbaqusINPInputT::ReadCoordinates (dArray2DT& coords)
{
	fAbaqus.Coordinates (coords);
}

void AbaqusINPInputT::ReadCoordinates (dArray2DT& coords, iArrayT& node_id)
{
  ReadCoordinates(coords);
  ReadNodeID(node_id);
}

int AbaqusINPInputT::NumGlobalElements (void) const
{
  cout << "\n\n AbaqusINPInputT::Not programmed to determine number of global elements \n\n";
	return -1;
}

int AbaqusINPInputT::NumElements (const StringT& name)
{
  StringT t = name;
  t.ToUpper();
  return fAbaqus.ElSetLength (t);
}

int AbaqusINPInputT::NumElementNodes (const StringT& name)
{
  StringT t = name;
  t.ToUpper ();
  return fAbaqus.NumElementNodesforSet (t);  
}

void AbaqusINPInputT::ReadAllElementMap (iArrayT& elemmap)
{
  cout << "\n\n AbaqusINPInputT::Not programmed to read all element map\n\n";
  elemmap = -1;
}

void AbaqusINPInputT::ReadGlobalElementMap (const StringT& name, iArrayT& elemmap)
{
  // element ids
  StringT t = name;
  t.ToUpper();
  fAbaqus.ElSet (t, elemmap);
}

void AbaqusINPInputT::ReadGlobalElementSet (const StringT& name, iArrayT& set)
{
	// must convert to index numbering
  cout << "\n\n AbaqusINPInputT::Not programmed to read element set\n\n";
  set = -1;
}

void AbaqusINPInputT::ReadConnectivity (const StringT& name, iArray2DT& connects)
{
  iArray2DT temp (connects.MajorDim(), connects.MinorDim());
  StringT t = name;
  t.ToUpper();
  if (!fAbaqus.Connectivity (t, temp))
  {
	  cout << "\n\n AbaqusINPInputT::Unable to read connectivity " << t << "\n\n";
	  throw ExceptionT::kDatabaseFail;
  }
  
  // convert node ids tags to index positions
  iArrayT node_ids (NumNodes());
  ReadNodeID (node_ids);
  int *p = temp.Pointer();
  int *c = connects.Pointer();
  for (int i=0; i < connects.Length(); i++)
  {
  	int index;
  	if (node_ids.HasValue (*p++, index))
  		*c++ = index;
  	else
  		throw ExceptionT::kOutOfRange;
  }
}

void AbaqusINPInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& code)
{
	code = fAbaqus.GeometryCode (name);
}


int AbaqusINPInputT::NumNodesInSet (const StringT& name)
{
  StringT t = name;
  t.ToUpper ();
  return fAbaqus.NodeSetLength (t);
}

void AbaqusINPInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
  iArrayT temp (nodes.Length());
  StringT t = name;
  t.ToUpper();
  fAbaqus.NodeSet (t, temp);

  // convert node ids tags to index positions
  iArrayT node_ids (NumNodes());
  ReadNodeID (node_ids);
  int *p = temp.Pointer();
  int *c = nodes.Pointer();
  for (int i=0; i < nodes.Length(); i++)
  {
  	int index;
  	if (node_ids.HasValue (*p++, index))
  		*c++ = index;
  	else
  		throw ExceptionT::kOutOfRange;
  }
}

int AbaqusINPInputT::NumSidesInSet (const StringT& anme) const
{
#pragma unused(anme)

  cout << "\n\n AbaqusINPInputT::Not programmed to read side sets\n\n";
  return 0;
}

StringT AbaqusINPInputT::SideSetGroupName (const StringT& name) const
{
#pragma unused(name)

  cout << "\n\n AbaqusINPInputT::Not programmed to read side sets\n\n";
  StringT elname ("");
  return elname; 
}

void AbaqusINPInputT::ReadSideSetLocal (const StringT& name, iArray2DT& sides) const
{
#pragma unused(name)
#pragma unused(sides)

  cout << "\n\n AbaqusINPInputT::Not programmed to read side sets\n\n";
  throw ExceptionT::kDatabaseFail;
}

void AbaqusINPInputT::ReadSideSetGlobal (const StringT& name, iArray2DT& sides) const
{
#pragma unused(name)
#pragma unused(sides)

  cout << "\n\n AbaqusINPInputT::Not programmed to read side sets\n\n";
  throw ExceptionT::kDatabaseFail;
}

/**************** PRIVATE *******************/

