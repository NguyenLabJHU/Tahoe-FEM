/* $Id: ModelManagerT.cpp,v 1.11 2002-01-09 12:18:35 paklein Exp $ */
/* created: sawimme July 2001 */

#include "ModelManagerT.h"
#include <ctype.h>

#include "ifstreamT.h"
#include "nVariArray2DT.h"

#include "TahoeInputT.h"
#include "ExodusInputT.h"
#include "PatranInputT.h"
#include "EnSightInputT.h"
#include "AbaqusInputT.h"
#include "InputFEASCIIT.h"

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

bool ModelManagerT::Initialize (ifstreamT& in, bool readonly)
{
	/* clear any existing parameters */
	Clear();

  StringT database;
  in >> fFormat;

  if (fFormat != IOBaseT::kTahoe)
    {
      in >> database;

      /* prepend full path name to database file */
      database.ToNativePathName();
      
      /* patch from input file */
      StringT path;
      path.FilePath (in.filename());

      /* prepend path */
      database.Prepend (path);
    }
  else
    database = "\0";

	fInputName = database;
 	if (!readonly)
 		return ScanModel(fInputName);
 	else
 		return true;
}

bool ModelManagerT::Initialize (const IOBaseT::FileTypeT format, const StringT& database)
{
	/* clear any existing parameters */
	Clear();

	fFormat = format;
	fInputName = database;
	return ScanModel(fInputName);
}

bool ModelManagerT::Initialize (void)
{
	/* clear any existing parameters */
	Clear();

  IOBaseT temp (cout);
  temp.InputFormats (cout);
  StringT database;
  cout << "\n Enter the Model Format Type: ";
  cin >> fFormat;
  if (fFormat != IOBaseT::kTahoe)
    {
      cout << "\n Enter the Model File Name: ";
      cin >> database;
    }
  else
    database = "\0";
	fInputName = database;
	return ScanModel(fInputName);
}

void ModelManagerT::EchoData (ostream& o) const
{
  IOBaseT temp (o);
  o << " Input format. . . . . . . . . . . . . . . . . . = " << fFormat  << '\n';
  temp.InputFormats (o);
  if (fFormat != IOBaseT::kTahoe)
    o << " Geometry file . . . . . . . . . . . . . . . . . = " << fInputName  << '\n';
}

bool ModelManagerT::RegisterNodes (int length, int dof)
{
  fCoordinateDimensions[0] = length;
  fCoordinateDimensions[1] = dof;
  return true;
}

bool ModelManagerT::RegisterElementGroup (const StringT& name, int numelems, int numelemnodes, GeometryT::CodeT code)
{
  if (!CheckName (fElementNames, name, "Element Group")) return false;
  
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
  if (!CheckName (fElementNames, name, "Node Set")) return false;
  
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
  if (!CheckName (fElementNames, name, "Side Set")) return false;
  
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
  if (!RegisterNodes (coords.MajorDim(), coords.MinorDim())) 
    return false;
  fCoordinates = coords;
  return true;
}

bool ModelManagerT::RegisterElementGroup (const StringT& name, iArray2DT& conn, GeometryT::CodeT code)
{
  if (!CheckName (fElementNames, name, "Element Group")) return false;
  
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
  if (!CheckName (fElementNames, name, "Node Set")) return false;
  
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
  if (!CheckName (fElementNames, name, "Side Set")) return false;
  
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

/* reads dimensions and numbered array, then offsets array */
bool ModelManagerT::RegisterNodes (ifstreamT& in)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterNodes(ifstreamT): count not open file");

  in2 >> fCoordinateDimensions[0] >> fCoordinateDimensions[1];
  fCoordinates.Allocate (fCoordinateDimensions[0], fCoordinateDimensions[1]);
  fCoordinates.ReadNumbered (in2);
  return true;
}

/* read dimensions and numbered array, then offsets array */
bool ModelManagerT::RegisterElementGroup (ifstreamT& in, const StringT& name, GeometryT::CodeT code)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterElementGroup(ifstreamT): count not open file");

  int length, numnodes;
  in2 >> length >> numnodes;
  iArray2DT temp (length, numnodes);
  temp.ReadNumbered (in2);
  temp += -1;
  return RegisterElementGroup (name, temp, code);
}

/* read dimensions and array, then offsets array */
bool ModelManagerT::RegisterNodeSet (ifstreamT& in, const StringT& name)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterNodeSet(ifstreamT): count not open file");

  int length;
  in2 >> length;
  if (length > 0)
    {
      iArrayT n (length);
      in2 >> n;
      n--;
      return RegisterNodeSet (name, n);
    }
  else
    return false;
}

/* read dimensions and array, then offsets array */
bool ModelManagerT::RegisterSideSet (ifstreamT& in, const StringT& name, bool local, int groupindex)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterSideSet (ifstreamT): count not open file");

  int length;
  in2 >> length;
  if (length > 0)
    {
      iArray2DT s (length, 2);
      in2 >> s;
      s += -1;
      return RegisterSideSet (name, s, local, groupindex);
    }
  else
    return false;
}

void ModelManagerT::ReadInlineCoordinates (ifstreamT& in)
{
  if (fFormat == IOBaseT::kTahoe)
    RegisterNodes (in);
}

void ModelManagerT::ElementBlockList (ifstreamT& in, iArrayT& indexes, iArrayT& matnums)
{
  /* number of blocks in element data */
  int num_blocks = 0;
  in >> num_blocks;
  fMessage << " Number of connectivity data blocks. . . . . . . = " << num_blocks << '\n';
  if (num_blocks < 1) throw eBadInputValue;

  indexes.Allocate (num_blocks);
  matnums.Allocate (num_blocks);
  for (int i=0; i < num_blocks; i++)
    {
      /* read mat id */
      in >> matnums[i];
      fMessage << "                   material number: " << matnums[i] << '\n';

      /* read element group name */
      StringT name;
      if (fFormat == IOBaseT::kTahoe)
	{
	  name = "ElementGroup";
	  name.Append (fNumElementSets+1);
	  RegisterElementGroup (in, name, GeometryT::kNone);
	}
      else
	{
	  in >> name;
	}
      fMessage << "                element block name: " << name << endl;

      /* check with model manager */
      indexes[i] = ElementGroupIndex (name);
      if (indexes[i] < 0)
	{
	  cout << "\n ModelManagerT::ElementBlockList: "
		   << " block " << i+1 << ": " << name << '\n';
	  cout << "      does not match Model database " << endl;
	  throw eBadInputValue;
	}
    }
}

void ModelManagerT::NodeSetList (ifstreamT& in, iArrayT& indexes)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      StringT name = "InlineNS";
      name.Append (fNumNodeSets+1);
      RegisterNodeSet (in, name);

      indexes.Allocate (1);
      indexes[0] = NodeSetIndex (name);

      /* account for no sets or all nodes */
      if (indexes[0] > -1)
	{
	  fMessage << " Node Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Node Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << indexes[0] << '\n';
	  fMessage << " Node Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fNodeSetDimensions[indexes[0]] << '\n';
	}
    }
  else
    {
      int num_sets;
      in >> num_sets;

      indexes.Allocate (num_sets);
      for (int i=0; i < num_sets; i++)
	{
	  StringT name;
	  in >> name;
	  indexes[i] = NodeSetIndex (name);

	  fMessage << " Node Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Node Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << indexes[0] << '\n';
	  fMessage << " Node Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fNodeSetDimensions[indexes[i]] << '\n';
	}
    }
}

void ModelManagerT::SideSetList (ifstreamT& in, iArrayT& indexes, bool multidatabasesets)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      bool local = true;

      int blockID;
      in >> blockID;
      int groupindex = ElementGroupIndex (blockID);

      StringT name = "InlineSS";
      name.Append (fNumSideSets+1);
      RegisterSideSet (in, name, local, groupindex);

      indexes.Allocate (1);
      indexes[0] = SideSetIndex (name);

      /* account for no sets */
      if (indexes[0] > -1)
	{
	  fMessage << " Side Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Side Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << indexes[0] << '\n';
	  fMessage << " Side Set Element Group Name . . . . . . . . . . = ";
	  fMessage << blockID << '\n';
	  fMessage << " Side Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fSideSetDimensions[indexes[0]] << '\n';
	}
    }
  else
    {
      int num_sets;
      if (multidatabasesets)
	in >> num_sets;
      else
	num_sets = 1;

      indexes.Allocate (num_sets);
      for (int i=0; i < num_sets; i++)
	{
	  StringT name;
	  in >> name;
	  indexes[i] = SideSetIndex (name);

	  fMessage << " Side Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Side Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << indexes[0] << '\n';
	  fMessage << " Side Set Element Group Name . . . . . . . . . . = ";
	  fMessage << fInput->SideSetGroupName (name) << '\n';
	  fMessage << " Side Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fSideSetDimensions[indexes[i]] << '\n';
	}
    }
}

/* return the total number of nodes, read node lists, integer data and double values */
int ModelManagerT::ReadCards (ifstreamT& in, ostream& out, ArrayT<iArrayT>& nodes, iArray2DT& data, dArrayT& value)
{
  const int numc = value.Length();
  if (data.MajorDim() != numc ||
      nodes.Length () != numc ) throw eSizeMismatch;
  data = 0;
  value = 0;

  /* account for text file name instead of data */
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, out, true, "ModelManagerT::ReadCards: could not open file");

  int count = 0;
  int *pd = data.Pointer();
  StringT ID;
  for (int i=0; i < numc; i++)
    {
      /* read node set name or node number */
      in2 >> ID;

      /* read nodes in set or create a set from node number */
      if (fFormat == IOBaseT::kTahoe)
	{
	  nodes[i].Allocate (1);
	  nodes[i] = atoi (ID) - 1; // offset
	  count++;
	}
      else
	{
	  int index = NodeSetIndex (ID);
	  if (index < 0) 
	    {
	      cout << "ModelManagerT::ReadCards, cannot find node set\n";
	      cout << "   asking for : " << ID << ".\n";
	      cout << "   available sets: \n";
	      for (int ig=0; ig < fNumNodeSets; ig++)
		    cout << "       " << fNodeSetNames[ig] << ".\n";
	      cout << endl;
	      throw eBadInputValue;
	    }
	  nodes[i] = NodeSet (index);
	  if (i == 0)
	    out << " Number of node sets . . . . . . . . . . . . . . = " 
		<< numc << "\n\n";
	  out << " Node Set Name . . . . . . . . . . . . . . . . . = ";
	  out << ID << '\n';
	  out << " Number of cards . . . . . . . . . . . . . . . . = ";
	  out << nodes[i].Length() << endl;
	  count += nodes[i].Length();
	}

      /* read card data */
      for (int j=0; j < data.MinorDim(); j++)
	in2 >> *pd++;
      in2 >> value[i];
    }  
  return count;
}

void ModelManagerT::ReadNumTractionLines (ifstreamT& in, int& numlines, int& numsets)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      in >> numlines;
      fMessage << " Number of traction BC's . . . . . . . . . . . . = " << numlines << '\n';
      
      if (numlines > 0)
	in >> numsets;
    }
  else
    {
      in >> numsets;
      numlines = numsets;
      fMessage << " Number of traction BC side sets . . . . . . . . = " << numsets << "\n\n";
    }
}

void ModelManagerT::ReadTractionSetData (ifstreamT& in, int& blockindex, int& setsize)
{
  if (fFormat == IOBaseT::kTahoe)
    in >> blockindex >> setsize;
  else
    {
      setsize = 1;
      // blockindex set later
    }
}

void ModelManagerT::ReadTractionSideSet (ifstreamT& in, int& blockindex, iArray2DT& localsides)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      localsides.Allocate (2,1);
      in >> localsides[0] >> localsides[1];
      // blockindex is already set
    }
  else
    {
      StringT name;
      int index;
      in >> name;
      index = SideSetIndex (name);
      if (index < 0) 
	{
	  cout << "\n ModelManagerT::ReadTractionSideSet: ";
	  cout << "Side Set Name is not found in registered list: ";
	  cout << name << "\n";
	  throw eBadInputValue;
	}
      blockindex = SideSetGroupIndex (index);
 
      /* read set */
      iArray2DT temp = SideSet (index);
      if (!IsSideSetLocal(index))
	SideSetGlobalToLocal (blockindex, localsides, temp);
      else
	localsides = temp;

      fMessage << " Database side set name. . . . . . . . . . . . . = ";
      fMessage << name << '\n';
      fMessage << " Number of traction BC cards . . . . . . . . . . = ";
      fMessage << localsides.MajorDim() << endl;
    }
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
	  cout << "\n ModelManagerT::Coordinates, coords not registered yet" << endl;
	  throw eGeneralFail;
	}
      else
	return; // do nothing, already loaded
    }
  fCoordinates.Allocate (fCoordinateDimensions[0], fCoordinateDimensions[1]); 
	if (!fInput) {
    	cout << "\n ModelManagerT::ReadCoordinates: input source is not initialized" << endl;
		throw eDatabaseFail;
	}
	fInput->ReadCoordinates (fCoordinates);
}

/* used to reduce 3D database coordinates (Patran, Abaqus, etc.) */
bool ModelManagerT::AreElements2D (void) const
{
  if (fCoordinateDimensions[1] < 3) return true;

  /* look over registered element sets */
  for (int i=0; i < fNumElementSets; i++)
    if (fElementCodes[i] == GeometryT::kPoint ||
	fElementCodes[i] == GeometryT::kTriangle ||
	fElementCodes[i] == GeometryT::kQuadrilateral )
      return true;

  return false;
}

void ModelManagerT::ElementGroupNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < names.Length(); i++)
    names[i] = fElementNames[i];
}

int ModelManagerT::ElementGroupIndex (const StringT& name) const
{
  // account for space padding at end of name
  int length1 = name.Length();
  for (int i=0; i < fNumElementSets; i++)
    {
      int length2 = fElementNames[i].Length();
      int length = (length1 < length2) ? length1 : length2;
      if (strncmp (name.Pointer(), fElementNames[i].Pointer(), length-1) == 0)
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
  ReadConnectivity (index);
  return fElementSets [index];
}

void ModelManagerT::ReadConnectivity (int index)
{
  if (index < 0 && index >= fNumElementSets)
    throw eOutOfRange;
  if (fElementSets[index].Length() == 0)
    {
      if (fFormat == IOBaseT::kTahoe)
	{
	  cout << "\n ModelManagerT::ReadConnectivity, elems not registered yet" << endl;
	  throw eGeneralFail;
	}
      fElementSets[index].Allocate (fElementLengths[index], fElementNodes[index]);
     if (!fInput) {
    	cout << "\n ModelManagerT::ReadConnectivity: input source is not initialized" << endl;
		throw eDatabaseFail;
      }
      fInput->ReadConnectivity (fElementNames[index], fElementSets[index]);
    }
}

const iArray2DT* ModelManagerT::ElementGroupPointer (int index) const
{
  if (index < 0 && index >= fNumElementSets)
    throw eOutOfRange;
  return &fElementSets[index];
}

void ModelManagerT::AllNodeMap (iArrayT& map)
{
	if (!fInput) {
    	cout << "\n ModelManagerT::AllNodeMap: input source is not initialized" << endl;
		throw eDatabaseFail;
	}
	fInput->ReadNodeMap (map);
}

void ModelManagerT::AllElementMap (iArrayT& map)
{
	if (!fInput) {
    	cout << "\n ModelManagerT::AllElementMap: input source is not initialized" << endl;
		throw eDatabaseFail;
	}
	fInput->ReadAllElementMap (map);
}

void ModelManagerT::ElementMap (StringT& name, iArrayT& map)
{
	if (!fInput) {
    	cout << "\n ModelManagerT::ElementMap: input source is not initialized" << endl;
		throw eDatabaseFail;
	}
	fInput->ReadGlobalElementMap (name, map);
}

void ModelManagerT::NodeSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < names.Length(); i++)
    names[i] = fNodeSetNames[i];
}

int ModelManagerT::NodeSetIndex (const StringT& name) const
{
  // account for space padding at end of name
  int length1 = name.Length();
  for (int i=0; i < fNumNodeSets; i++)
    {
      int length2 = fNodeSetNames[i].Length();
      int length = (length1 < length2) ? length1 : length2;
      if (strncmp (name.Pointer(), fNodeSetNames[i].Pointer(), length-1) == 0)
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
	  cout << "\n ModelManagerT::NodeSet, set not registered yet" << endl;
	  throw eGeneralFail;
	}
      fNodeSets[index].Allocate (fNodeSetDimensions[index]);
      if (!fInput) {
    	cout << "\n ModelManagerT::NodeSet: input source is not initialized" << endl;
		throw eDatabaseFail;
	  }
      fInput->ReadNodeSet (fNodeSetNames[index], fNodeSets[index]);
    }
  return fNodeSets [index];
}

void ModelManagerT::ManyNodeSets (const iArrayT& indexes, iArrayT& nodes)
{
  iAutoArrayT temp;
  iArrayT tn;
  for (int i=0; i < indexes.Length(); i++)
    {
      tn = NodeSet (indexes[i]);
      temp.AppendUnique(tn);
    }

  nodes.Allocate (temp.Length());
  nodes.CopyPart (0, temp, 0, temp.Length());
  nodes.SortAscending ();
}

void ModelManagerT::SideSetNames (ArrayT<StringT>& names) const
{
  for (int i=0; i < names.Length(); i++)
    names[i] = fSideSetNames[i];
}

int ModelManagerT::SideSetIndex (const StringT& name) const
{
  // account for space padding at end of name
  int length1 = name.Length();
  for (int i=0; i < fNumSideSets; i++)
    {
      int length2 = fSideSetNames[i].Length();
      int length = (length1 < length2) ? length1 : length2;
      if (strncmp (name.Pointer(), fSideSetNames[i].Pointer(), length-1) == 0)
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
	  cout << "\n ModelManagerT::SideSet, set not registered yet" << endl;
	  throw eGeneralFail;
	}
      fSideSets[index].Allocate (fSideSetDimensions[index], 2);
      if (!fInput) {
    	cout << "\n ModelManagerT::SideSet: input source is not initialized" << endl;
		throw eDatabaseFail;
	  }
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
  cout << "\n ModelManagerT::SideSetGlobalToLocal not implemented" << endl;
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
  int newnodes = newcoords.MajorDim();
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

void ModelManagerT::AdjustCoordinatesto2D (void)
{
  if (fCoordinateDimensions[1] != 3) return;

  /* make sure coordinates are already loaded */
  ReadCoordinates ();

  /* copy first two dimensions */
  dArray2DT temp (fCoordinateDimensions[0], 2);
  double *pt = temp.Pointer();
  double *pc = fCoordinates.Pointer();
  for (int i=0; i < fCoordinateDimensions[0]; i++)
    {
      for (int j=0; j < 2; j++)
	*pt++ = *pc++;
      *pc++;
    }
  
  /* overwrite registered values */
  RegisterNodes (temp);
}

bool ModelManagerT::RegisterVariElements (const StringT& name, nVariArray2DT<int>& conn, 
					  GeometryT::CodeT code, int numelemnodes,
					  int headroom)
{
  if (!CheckName (fElementNames, name, "Element Group")) return false;
  
  fElementNames.Append (name);
  fElementLengths.Append (0);
  fElementNodes.Append (numelemnodes);
  fElementCodes.Append (code);

  int index = fElementSets.Length();
  iArray2DT temp;
  fElementSets.Append (temp);

  fNumElementSets++;
  if (fNumElementSets != fElementNames.Length() ||
      fNumElementSets != fElementLengths.Length() ||
      fNumElementSets != fElementNodes.Length() ||
      fNumElementSets != fElementCodes.Length() ||
      fNumElementSets != fElementSets.Length() )
    return false;

  conn.SetWard (headroom, fElementSets [index], numelemnodes);

  return true;
}

/* call this function after the connectivity has been changed by outside classes */
void ModelManagerT::UpdateConnectivity (int index, const iArray2DT& connects)
{
  if (index < 0 && index >= fNumElementSets) throw eOutOfRange;

  fElementSets[index] = connects;

  fElementLengths[index] = fElementSets[index].MajorDim();
  fElementNodes[index] = fElementSets[index].MinorDim();
}

void ModelManagerT::AddElement (int index, const iArray2DT& connects, iArrayT& new_elem_tags, int& numelems)
{
  if (index < 0 && index >= fNumElementSets) throw eOutOfRange;
  if (connects.MajorDim() != new_elem_tags.Length() ||
      connects.MinorDim() != fElementSets[index].MinorDim() ) throw eSizeMismatch;

  /* set new elem tags to the old last elem number + 1 */
  new_elem_tags.SetValueToPosition ();
  new_elem_tags += fElementSets[index].MajorDim();

  /* reset the number of elements */
  int newelems = connects.MajorDim();
  fElementLengths[index] += newelems;
  numelems = fElementLengths[index];

  /* reallocate */
  fElementSets[index].Resize (newelems);

  /* copy in */
  int *pc = connects.Pointer();
  for (int i=0; i < newelems; i++, pc += connects.MinorDim())
    fElementSets[index].SetRow (new_elem_tags[i], pc);
}

void ModelManagerT::CloseModel (void)
{
  if (fInput) fInput->Close ();
  delete fInput;
  fInput = NULL;
}

ifstreamT& ModelManagerT::OpenExternal (ifstreamT& in, ifstreamT& in2, ostream& out, bool verbose, const char* fail) const
{
  /* external files are only allowed with inline text */
  if (fFormat != IOBaseT::kTahoe)
    return in;

	/* check for external file */
	char nextchar = in.next_char();
	if (isdigit(nextchar))
		return in;
	else
	{
		/* open external file */
		StringT file;
		in >> file;
		if (verbose) out << " external file: " << file << '\n';
		file.ToNativePathName();

		/* path to source file */
		StringT path;
		path.FilePath(in.filename());
		file.Prepend(path);
			
		/* open stream */
		in2.open(file);
		if (!in2.is_open())
		{
			if (verbose && fail) cout << "\n " << fail << ": " << file << endl;
			throw eBadInputValue;
		}

		/* set comments */
		if (in.skip_comments()) in2.set_marker(in.comment_marker());

		return in2;
	}  
}

/*********** PRIVATE **************/

bool ModelManagerT::ScanModel (const StringT& database)
{
  switch (fFormat)
    {
    case IOBaseT::kTahoe:
      /* do nothing, arrays will be registered via ElementBaseT and NodeManager */
      fInput = NULL;
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
    case IOBaseT::kTahoeResults:
      fInput = new InputFEASCIIT (fMessage);
      break;
    default:
      {
	cout << "\n ModelManagerT::Unsupported model format: " 
		 << fFormat << endl;
	throw eGeneralFail;
      }
    }

	if (fFormat != IOBaseT::kTahoe)
	{
		if (!fInput) 
		{
			fMessage << "\n ModelManagerT::Scan Model fInput not set" << endl;
			return false;
		}

		fInput->Close();
		if (!fInput->Open(database))
		{
			fMessage << "\n ModelManagerT::ScanModel: error opening database file \"" 
			         << database << '\"' << endl;
			return false;
		}
      
		fCoordinateDimensions [0] = fInput->NumNodes ();
		fCoordinateDimensions [1] = fInput->NumDimensions ();

		if (!ScanElements ())
		{
			fMessage << "\n ModelManagerT::ScanModel: Error Registering Elements" << endl;
			return false;
		}
 
		if (!ScanNodeSets ())
		{
			fMessage << "\n ModelManagerT::ScanModel: Error Registering NodeSets" << endl;
			return false;
		}

		if (!ScanSideSets ())
		{
			fMessage << "\n ModelManagerT::ScanModel: Error Registering SideSets" << endl;
			return false;
		}
	}
	
	/* success */
	return true;
}

bool ModelManagerT::ScanElements (void)
{
	/* check if input has been initialized */
	if (!fInput) return false;
	
  fNumElementSets = fInput->NumElementGroups ();
  fElementLengths.Allocate (fNumElementSets);
  fElementNodes.Allocate (fNumElementSets);
  fElementNames.Allocate (fNumElementSets);
  fElementCodes.Allocate (fNumElementSets);
  fElementSets.Allocate (fNumElementSets);

  if (fNumElementSets > 0)
    {
      fInput->ElementGroupNames (fElementNames);
      for (int e=0; e < fNumElementSets; e++)
	{
	  fElementLengths[e] = fInput->NumElements (fElementNames[e]);
	  fElementNodes[e] = fInput->NumElementNodes (fElementNames[e]);
	  fInput->ReadGeometryCode (fElementNames[e], fElementCodes[e]);
	}
    }
  return true;
}

bool ModelManagerT::ScanNodeSets (void)
{
	/* check if input has been initialized */
	if (!fInput) return false;

  fNumNodeSets = fInput->NumNodeSets();
  fNodeSetNames.Allocate (fNumNodeSets);
  fNodeSetDimensions.Allocate (fNumNodeSets);
  fNodeSets.Allocate (fNumNodeSets);
  if (fNumNodeSets > 0)
    {
      fInput->NodeSetNames (fNodeSetNames);
      for (int i=0; i < fNumNodeSets; i++)
	fNodeSetDimensions[i] = fInput->NumNodesInSet (fNodeSetNames[i]);
    }
  return true;
}

bool ModelManagerT::ScanSideSets (void)
{
	/* check if input has been initialized */
	if (!fInput) return false;

  fNumSideSets = fInput->NumSideSets();
  fSideSetNames.Allocate (fNumSideSets);
  fSideSetDimensions.Allocate (fNumSideSets);
  fSideSets.Allocate (fNumSideSets);
  fSideSetIsLocal.Allocate (fNumSideSets);
  fSideSetGroupIndex.Allocate (fNumSideSets);

  if (fNumSideSets > 0)
    {
      fInput->SideSetNames (fSideSetNames);
      bool t = fInput->AreSideSetsLocal ();
      fSideSetIsLocal = t;
      fSideSetGroupIndex = -1;

      for (int i=0; i < fNumSideSets; i++)
	{
	  fSideSetDimensions[i] = fInput->NumSidesInSet (fSideSetNames[i]);
	  if (fSideSetIsLocal[i])
	    {
	      StringT name = fInput->SideSetGroupName (fSideSetNames[i]);
	      fSideSetGroupIndex[i] = ElementGroupIndex (name);
	    }
	}
    }
  return true;
}


bool ModelManagerT::CheckName (const ArrayT<StringT>& list, const StringT& name, const char *type) const
{
  // account for space padding at end of name
  int l1 = name.Length();

  for (int i=0; i < list.Length(); i++)
    {
      int l2 = list[i].Length();
      int l = (l1 < l2) ? l1 : l2;
      if (strncmp (list[i].Pointer(), name.Pointer(), l-1) == 0)
	{
	  fMessage << "\nModelManagerT::CheckName\n";
	  fMessage << "   " << type << " already has a registered set called " << name << "\n\n";
	  fMessage << "  Sets: \n";
	  for (int j=0; j < list.Length(); j++)
	    fMessage << "       " << list[i] << "\n";
	  fMessage << "\n";
	  return false;
	}
    }
  return true;
}

/* clear database parameters */
void ModelManagerT::Clear(void)
{
	/* close the input */
	delete fInput;
	fInput = NULL;

	//TO DO: clear all memory and set arrays back to empty
}
