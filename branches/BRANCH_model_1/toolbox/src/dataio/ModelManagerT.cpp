/* $Id: ModelManagerT.cpp,v 1.4.2.3 2001-10-04 22:06:28 sawimme Exp $ */
/* created: sawimme July 2001 */

#include "ModelManagerT.h"
#include "ifstreamT.h"

#include "TahoeInputT.h"
#include "ExodusInputT.h"
#include "PatranInputT.h"
#include "EnSightInputT.h"
#include "AbaqusInputT.h"

ModelManagerT::ModelManagerT (ostream& message) :
  fCoordinateDimensions (2),
  fNumElementSets (0),
  fNumNodeSets (0),
  fNumSideSets (0),
  fInput (NULL),
  fInputName (""),
  fMessage (message),
  fFormat (IOBaseT::kTahoe)
{
  fCoordinateDimensions = -1;
}

ModelManagerT::~ModelManagerT (void)
{
  delete fInput;
}

void ModelManagerT::Initialize (ifstreamT& in)
{
  StringT database;
  in >> fFormat;
  in >> database;
  ScanModel (database);
}

void ModelManagerT::Initialize (const IOBaseT::FileTypeT format, const StringT& database)
{
  fFormat = format;
  ScanModel (database);
}

void ModelManagerT::Initialize (void)
{
  IOBaseT temp (cout);
  temp.InputFormats (cout);
  StringT database;
  cout << "\n Enter the Model Format Type: ";
  cin >> fFormat;
  cout << "\n Enter the Model File Name: ";
  cin >> database;
  ScanModel (database);
}

bool ModelManagerT::RegisterNodes (int length, int dof)
{
  fCoordinateDimensions[0] = length;
  fCoordinateDimensions[1] = dof;
  return true;
}

bool ModelManagerT::RegisterElementGroup (const StringT& name, int numelems, int numelemnodes, GeometryT::CodeT code)
{
  fElementNames.Append (name);
  fElementLengths.Append (numelems);
  fElementNodes.Append (numelemnodes);
  fElementCodes.Append (code);

  iArray2DT temp;
  fElementSets.Append (temp);

  fNumElementSets++;
  if (fNumElementSets != fElementNames.Length() ||
      fNumElementSets != fElementLengths.Length() ||
      fNumElementSets != fElementNodes.Length() ||
      fNumElementSets != fElementCodes.Length() ||
      fNumElementSets != fElementSets.Length() )
    return false;
  return true;
}

bool ModelManagerT::RegisterNodeSet (const StringT& name, int length)
{
  fNodeSetNames.Append (name);
  fNodeSetDimensions.Append (length);

  iArrayT temp;
  fNodeSets.Append (temp);

  fNumNodeSets++;
  if (fNumNodeSets != fNodeSetNames.Length() ||
      fNumNodeSets != fNodeSetDimensions.Length() ||
      fNumNodeSets != fNodeSets.Length() )
    return false;
  return true;
}

bool ModelManagerT::RegisterSideSet (const StringT& name, int length, bool local, int elemgroupindex)
{
  fSideSetNames.Append (name);
  fSideSetDimensions.Append (length);
  fSideSetIsLocal.Append (local);
  if (local)
    fSideSetGroupIndex.Append (elemgroupindex);
  else
    fSideSetGroupIndex.Append (-1);

  iArray2DT temp;
  fSideSets.Append (temp);

  fNumSideSets++;
  if (fNumSideSets != fSideSetNames.Length() ||
      fNumSideSets != fSideSetDimensions.Length() ||
      fNumSideSets != fSideSets.Length() ||
      fNumSideSets != fSideSetIsLocal.Length() ||
      fNumSideSets != fSideSetGroupIndex.Length() )
    return false;
  return true;
}

bool ModelManagerT::RegisterNodes (dArray2DT& coords)
{
  fCoordinateDimensions[0] = coords.MajorDim();
  fCoordinateDimensions[1] = coords.MinorDim();
  fCoordinates = coords;
  return true;
}

bool ModelManagerT::RegisterElementGroup (const StringT& name, iArray2DT& conn, GeometryT::CodeT code)
{
  fElementNames.Append (name);
  fElementLengths.Append (conn.MajorDim());
  fElementNodes.Append (conn.MinorDim());
  fElementCodes.Append (code);
  fElementSets.Append (conn);

  fNumElementSets++;
  if (fNumElementSets != fElementNames.Length() ||
      fNumElementSets != fElementLengths.Length() ||
      fNumElementSets != fElementNodes.Length() ||
      fNumElementSets != fElementCodes.Length() ||
      fNumElementSets != fElementSets.Length() )
    return false;
  return true;
}

bool ModelManagerT::RegisterNodeSet (const StringT& name, iArrayT& set)
{
  fNodeSetNames.Append (name);
  fNodeSetDimensions.Append (set.Length());
  fNodeSets.Append (set);

  fNumNodeSets++;
  if (fNumNodeSets != fNodeSetNames.Length() ||
      fNumNodeSets != fNodeSetDimensions.Length() ||
      fNumNodeSets != fNodeSets.Length() )
    return false;
  return true;
}

bool ModelManagerT::RegisterSideSet (const StringT& name, iArray2DT& set, bool local, int groupindex)
{
  fSideSetNames.Append (name);
  fSideSetDimensions.Append (set.MajorDim());
  fSideSets.Append (set);
  fSideSetIsLocal.Append (local);
  if (local)
    fSideSetGroupIndex.Append (groupindex);
  else
    fSideSetGroupIndex.Append (-1);

  fNumSideSets++;
  if (fNumSideSets != fSideSetNames.Length() ||
      fNumSideSets != fSideSetDimensions.Length() ||
      fNumSideSets != fSideSets.Length() ||
      fNumSideSets != fSideSetIsLocal.Length() ||
      fNumSideSets != fSideSetGroupIndex.Length() )
    return false;
  return true;
}

void ModelManagerT::CoordinateDimensions (int& length, int& dof) const
{
  length = fCoordinateDimensions[0];
  dof = fCoordinateDimensions[1];
}

const dArray2DT& ModelManagerT::CoordinateReference (void) const
{
  return fCoordinates;
}

const dArray2DT& ModelManagerT::Coordinates (void)
{
  if (fCoordinates.Length() == 0)
    ReadCoordinates ();
  return fCoordinates;
}

void ModelManagerT::ReadCoordinates (void)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      if (fCoordinates.Length() == 0)
	{
	  fMessage << "\n\nModelManagerT::Coordinates, coords not registered yet\n\n";
	  throw eGeneralFail;
	}
      else
	return; // do nothing, already loaded
    }
  fCoordinates.Allocate (fCoordinateDimensions[0], fCoordinateDimensions[1]); 
  if (!fInput) throw eGeneralFail;
  fInput->ReadCoordinates (fCoordinates);
}

void ModelManagerT::ElementGroupNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < names.Length(); i++)
    names[i] = fElementNames[i];
}

int ModelManagerT::ElementGroupIndex (const StringT& name) const
{
  for (int i=0; i < fNumElementSets; i++)
    {
      int length = fElementNames[i].Length();
      if (strncmp (name.Pointer(), fElementNames[i].Pointer(), length) == 0)
	return i;
    }
  return -1;
}

void ModelManagerT::ElementGroupDimensions (int index, int& numelems, int& numelemnodes) const
{
  numelems = -1;
  numelemnodes = -1;
  if (index < 0 && index >= fNumElementSets)
    return;
  numelems = fElementLengths[index];
  numelemnodes = fElementNodes[index];
}

GeometryT::CodeT ModelManagerT::ElementGroupGeometry (int index) const
{
  if (index < 0 && index >= fNumElementSets)
    return GeometryT::kNone;
  return fElementCodes[index];
}

const iArray2DT& ModelManagerT::ElementGroup (int index)
{
  if (index < 0 && index >= fNumElementSets)
    throw eOutOfRange;
  if (fElementSets[index].Length() == 0)
    {
      if (fFormat == IOBaseT::kTahoe)
	{
	  fMessage << "\n\nModelManagerT::ElementGroup, elems not registered yet\n\n";
	  throw eGeneralFail;
	}
      fElementSets[index].Allocate (fElementLengths[index], fElementNodes[index]);
      if (!fInput) throw eGeneralFail;
      fInput->ReadConnectivity (fElementNames[index], fElementSets[index]);
    }
  return fElementSets [index];
}

void ModelManagerT::AllNodeMap (iArrayT& map)
{
  if (!fInput) throw eGeneralFail;
  fInput->ReadNodeMap (map);
}

void ModelManagerT::AllElementMap (iArrayT& map)
{
  if (!fInput) throw eGeneralFail;
  fInput->ReadAllElementMap (map);
}

void ModelManagerT::ElementMap (StringT& name, iArrayT& map)
{
  if (!fInput) throw eGeneralFail;
  fInput->ReadGlobalElementMap (name, map);
}

void ModelManagerT::NodeSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < names.Length(); i++)
    names[i] = fNodeSetNames[i];
}

int ModelManagerT::NodeSetIndex (const StringT& name) const
{
  for (int i=0; i < fNumNodeSets; i++)
    {
      int length = fNodeSetNames[i].Length();
      if (strncmp (name.Pointer(), fNodeSetNames[i].Pointer(), length) == 0)
	return i;
    }
  return -1;
}

int ModelManagerT::NodeSetLength (int index) const
{
  if (index < 0 && index >= fNumNodeSets)
    return -1;
  return fNodeSetDimensions [index];
}

const iArrayT& ModelManagerT::NodeSet (int index)
{
  if (index < 0 && index >= fNumNodeSets)
    throw eOutOfRange;
  if (fNodeSets[index].Length() == 0)
    {
      if (fFormat == IOBaseT::kTahoe)
	{
	  fMessage << "\n\nModelManagerT::NodeSet, set not registered yet\n\n";
	  throw eGeneralFail;
	}
      fNodeSets[index].Allocate (fNodeSetDimensions[index]);
      if (!fInput) throw eGeneralFail;
      fInput->ReadNodeSet (fNodeSetNames[index], fNodeSets[index]);
    }
  return fNodeSets [index];
}

void ModelManagerT::SideSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < names.Length(); i++)
    names[i] = fSideSetNames[i];
}

int ModelManagerT::SideSetIndex (const StringT& name) const
{
  for (int i=0; i < fNumSideSets; i++)
    {
      int length = fSideSetNames[i].Length();
      if (strncmp (name.Pointer(), fSideSetNames[i].Pointer(), length) == 0)
	return i;
    }
  return -1;
}

int ModelManagerT::SideSetLength (int index) const
{
  if (index < 0 && index >= fNumSideSets)
    return -1;
  return fSideSetDimensions [index];
}

const iArray2DT& ModelManagerT::SideSet (int index) const
{
  if (index < 0 && index >= fNumSideSets)
    throw eOutOfRange;
  if (fSideSets[index].Length() == 0)
    {
      if (fFormat == IOBaseT::kTahoe)
	{
	  fMessage << "\n\nModelManagerT::SideSet, set not registered yet\n\n";
	  throw eGeneralFail;
	}
      fSideSets[index].Allocate (fSideSetDimensions[index], 2);
      if (!fInput) throw eGeneralFail;
      if (fSideSetIsLocal[index])
	fInput->ReadSideSetLocal (fSideSetNames[index], fSideSets[index]);
      else
	fInput->ReadSideSetGlobal (fSideSetNames[index], fSideSets[index]);
    }
  return fSideSets [index];
}

bool ModelManagerT::IsSideSetLocal (int index) const
{
  if (index < 0 && index >= fNumSideSets)
    throw eOutOfRange;
  return fSideSetIsLocal [index];
}

int ModelManagerT::SideSetGroupIndex (int sidesetindex) const
{
  if (sidesetindex < 0 && sidesetindex >= fNumSideSets)
    throw eOutOfRange;

  /* need to check if group index < 0, then have global side set and
     need to determine correct group index */

  return fSideSetGroupIndex [sidesetindex];
}

void ModelManagerT::SideSetLocalToGlobal (const int localelemindex, const iArray2DT& local, iArray2DT& global)
{
  int offset = 0;
  for (int i=0; i < localelemindex; i++)
    offset += fElementLengths[i];

  global = local;
  int *pelem = global.Pointer();
  for (int j=0; j < global.MajorDim(); j++, pelem += 2)
    *pelem += offset;
}

void ModelManagerT::SideSetGlobalToLocal (int& localelemindex, iArray2DT& local, const iArray2DT& global)
{
#pragma unused(localelemindex)
#pragma unused(local)
#pragma unused(global)
  cout << "\n\n ModelManagerT not programmed SideSetGlobalToLocal\n\n";
  throw eGeneralFail;
}

void ModelManagerT::AddNodes (const dArray2DT& newcoords, iArrayT& new_node_tags, int& numnodes)
{
  if (newcoords.MajorDim() != new_node_tags.Length() ||
      newcoords.MinorDim() != fCoordinates.MinorDim() ) throw eSizeMismatch;

  /* set new node tags to the old last node number + 1 */
  new_node_tags.SetValueToPosition ();
  new_node_tags += fCoordinateDimensions[0];

  /* reset the number of nodes */
  int newnodes = newcoords.MinorDim();
  fCoordinateDimensions[0] += newnodes;
  numnodes = fCoordinateDimensions[0];

  /* reallocate */
  fCoordinates.Resize (fCoordinateDimensions[0]);

  /* copy in */
  double *pc = newcoords.Pointer();
  for (int i=0; i < newnodes; i++, pc += newcoords.MinorDim())
    fCoordinates.SetRow (new_node_tags[i], pc);
}

void ModelManagerT::DuplicateNodes (const iArrayT& nodes, iArrayT& new_node_tags, int& numnodes)
{
  if (nodes.Length() != new_node_tags.Length()) throw eSizeMismatch;

  /* set new node tags to the old last node number + 1 */
  new_node_tags.SetValueToPosition ();
  new_node_tags += fCoordinateDimensions[0];

  /* reset the number of nodes */
  int newnodes = nodes.Length();
  fCoordinateDimensions[0] += newnodes;
  numnodes = fCoordinateDimensions[0];

  /* reallocate */
  fCoordinates.Resize (fCoordinateDimensions[0]);

  /* copy in */
  for (int i=0; i < nodes.Length(); i++)
    fCoordinates.CopyRowFromRow (new_node_tags[i], nodes[i]);
}

void ModelManagerT::CloseModel (void)
{
  if (fInput) fInput->Close ();
  delete fInput;
  fInput = NULL;
}

/*********** PRIVATE **************/

void ModelManagerT::ScanModel (const StringT& database)
{
  switch (fFormat)
    {
    case IOBaseT::kTahoe:
      /* do nothing, arrays will be registered via ElementBaseT and NodeManager */
      break;
      case IOBaseT::kTahoeII:
      fInput = new TahoeInputT (fMessage);
      break;
    case IOBaseT::kEnSight:
      fInput = new EnSightInputT (fMessage, false);
      break;
    case IOBaseT::kEnSightBinary:
      fInput = new EnSightInputT (fMessage, true);
      break;
    case IOBaseT::kExodusII:
      fInput = new ExodusInputT (fMessage);
      break;
    case IOBaseT::kPatranNeutral:
      fInput = new PatranInputT (fMessage);
      break;
    case IOBaseT::kAbaqus:
    case IOBaseT::kAbaqusBinary:
      fInput = new AbaqusInputT (fMessage);
      break;
    default:
      {
	fMessage << "\n\nModelManagerT::Unsupported model format. " 
		 << fFormat << "\n\n";
	throw eGeneralFail;
      }
    }

  if (fFormat != IOBaseT::kTahoe)
    {
      if (!fInput) throw eGeneralFail;
      fInput->Close ();
      fInput->Open (database);
      
      fCoordinateDimensions [0] = fInput->NumNodes ();
      fCoordinateDimensions [1] = fInput->NumDimensions ();

      if (!ScanElements ())
	{
	  fMessage << "\n\nModelManagerT::ScanModel: Error Registering Elements.\n\n";
	  throw eGeneralFail;
	}

      if (!ScanNodeSets ())
	{
	  fMessage << "\n\nModelManagerT::ScanModel: Error Registering NodeSets.\n\n";
	  throw eGeneralFail;
	}

      if (!ScanSideSets ())
	{
	  fMessage << "\n\nModelManagerT::ScanModel: Error Registering SideSets.\n\n";
	  throw eGeneralFail;
	}

      fInputName = database;
    }
}

bool ModelManagerT::ScanElements (void)
{
  fNumElementSets = fInput->NumElementGroups ();
  fElementLengths.Allocate (fNumElementSets);
  fElementNodes.Allocate (fNumElementSets);
  fElementNames.Allocate (fNumElementSets);
  fElementCodes.Allocate (fNumElementSets);
  fElementSets.Allocate (fNumElementSets);

  fInput->ElementGroupNames (fElementNames);
  for (int e=0; e < fNumElementSets; e++)
    {
      fElementLengths[e] = fInput->NumElements (fElementNames[e]);
      fElementNodes[e] = fInput->NumElementNodes (fElementNames[e]);
      fInput->ReadGeometryCode (fElementNames[e], fElementCodes[e]);
    }
  return true;
}

bool ModelManagerT::ScanNodeSets (void)
{
  fNumNodeSets = fInput->NumNodeSets();
  fNodeSetNames.Allocate (fNumNodeSets);
  fNodeSetDimensions.Allocate (fNumNodeSets);
  fNodeSets.Allocate (fNumNodeSets);
  fInput->NodeSetNames (fNodeSetNames);
  for (int i=0; i < fNumNodeSets; i++)
    fNodeSetDimensions[i] = fInput->NumNodesInSet (fNodeSetNames[i]);
  return true;
}

bool ModelManagerT::ScanSideSets (void)
{
  fNumSideSets = fInput->NumSideSets();
  fSideSetNames.Allocate (fNumSideSets);
  fSideSetDimensions.Allocate (fNumSideSets);
  fSideSets.Allocate (fNumSideSets);
  fSideSetIsLocal.Allocate (fNumSideSets);
  fSideSetGroupIndex.Allocate (fNumSideSets);

  fInput->SideSetNames (fSideSetNames);
  bool t = fInput->AreSideSetsLocal ();
  fSideSetIsLocal = t;
  fSideSetGroupIndex = -1;

  for (int i=0; i < fNumSideSets; i++)
    {
      fSideSetDimensions[i] = fInput->NumSidesInSet (fSideSetNames[i]);
      if (fSideSetIsLocal[i])
	fInput->SideSetGroupIndex (fSideSetNames[i]);
    }
  return true;
}
