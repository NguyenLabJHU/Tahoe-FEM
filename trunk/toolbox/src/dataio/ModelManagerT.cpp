/* $Id: ModelManagerT.cpp,v 1.15 2002-02-04 19:25:19 paklein Exp $ */
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

bool ModelManagerT::Initialize (const IOBaseT::FileTypeT format, const StringT& database, bool scan_model)
{
	/* clear any existing parameters */
	Clear();

	fFormat = format;
	fInputName = database;
	if (scan_model)
		return ScanModel(fInputName);
	else
		return true;
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
      database.ToNativePathName();
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

bool ModelManagerT::RegisterElementGroup (const StringT& ID, int numelems, 
	int numelemnodes, GeometryT::CodeT code)
{
  if (!CheckID (fElementNames, ID, "Element Group")) return false;
  
  fElementNames.Append (ID);
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

bool ModelManagerT::RegisterNodeSet (const StringT& ID, int length)
{
  if (!CheckID (fNodeSetNames, ID, "Node Set")) return false;
  
  fNodeSetNames.Append (ID);
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

bool ModelManagerT::RegisterSideSet (const StringT& ss_ID, int length, bool local, 
	const StringT& element_ID)
{
  if (!CheckID (fSideSetNames, ss_ID, "Side Set")) return false;
  
  fSideSetNames.Append (ss_ID);
  fSideSetDimensions.Append (length);
  fSideSetIsLocal.Append (local);
  if (local)
  {
  	int index = ElementGroupIndex(element_ID);
  	if (index == -1) {
  		cout << "\n ModelManagerT::RegisterSideSet: element ID not found: " << element_ID << endl;
  		throw eOutOfRange;
  	}
    fSideSetGroupIndex.Append (index);
  }
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

bool ModelManagerT::RegisterNodes (const dArray2DT& coords)
{
  if (!RegisterNodes (coords.MajorDim(), coords.MinorDim())) return false;
  fCoordinates = coords;
  return true;
}

bool ModelManagerT::RegisterElementGroup (const StringT& ID, const iArray2DT& conn, 
	GeometryT::CodeT code)
{
  if (!CheckID (fElementNames, ID, "Element Group")) return false;
  
  fElementNames.Append (ID);
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

bool ModelManagerT::RegisterNodeSet (const StringT& ID, const iArrayT& set)
{
  if (!CheckID (fNodeSetNames, ID, "Node Set")) return false;
  
  fNodeSetNames.Append (ID);
  fNodeSetDimensions.Append (set.Length());
  fNodeSets.Append (set);

  fNumNodeSets++;
  if (fNumNodeSets != fNodeSetNames.Length() ||
      fNumNodeSets != fNodeSetDimensions.Length() ||
      fNumNodeSets != fNodeSets.Length() )
    return false;
  return true;
}

bool ModelManagerT::RegisterSideSet (const StringT& ss_ID, const iArray2DT& set, 
	bool local, const StringT& element_ID)
{
	bool result = RegisterSideSet(ss_ID, set.MajorDim(), local, element_ID);
	fSideSets.Append (set);
	return result && fNumSideSets == fSideSets.Length();
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
bool ModelManagerT::RegisterElementGroup (ifstreamT& in, const StringT& ID, GeometryT::CodeT code)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterElementGroup(ifstreamT): count not open file");

  int length, numnodes;
  in2 >> length >> numnodes;
  iArray2DT temp (length, numnodes);
  temp.ReadNumbered (in2);
  temp += -1;
  return RegisterElementGroup (ID, temp, code);
}

/* read dimensions and array, then offsets array */
bool ModelManagerT::RegisterNodeSet (ifstreamT& in, const StringT& ID)
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
      return RegisterNodeSet (ID, n);
    }
  else
    return false;
}

/* read dimensions and array, then offsets array */
bool ModelManagerT::RegisterSideSet (ifstreamT& in, const StringT& ss_ID, bool local, 
	const StringT& element_ID)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterSideSet (ifstreamT): count not open file");

  int length;
  in2 >> length;
  if (length > 0)
    {
      iArray2DT s (length, 2);
      in2 >> s;
      s--;
      return RegisterSideSet (ss_ID, s, local, element_ID);
    }
  else
    return false;
}

void ModelManagerT::ReadInlineCoordinates (ifstreamT& in)
{
  if (fFormat == IOBaseT::kTahoe)
    RegisterNodes (in);
}

void ModelManagerT::ElementBlockList (ifstreamT& in, ArrayT<StringT>& ID, iArrayT& matnums)
{
  /* number of blocks in element data */
  int num_blocks = 0;
  in >> num_blocks;
  fMessage << " Number of connectivity data blocks. . . . . . . = " << num_blocks << '\n';
  if (num_blocks < 1) throw eBadInputValue;

  ID.Allocate (num_blocks);
  matnums.Allocate (num_blocks);
  for (int i=0; i < num_blocks; i++)
    {
      /* read mat id */
      in >> matnums[i];
      fMessage << "                   material number: " << matnums[i] << '\n';

      /* read element group name */
      StringT& name = ID[i];
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

    }
}

void ModelManagerT::NodeSetList (ifstreamT& in, ArrayT<StringT>& ID)
{
	if (fFormat == IOBaseT::kTahoe)
	{
		/* read set */
		StringT name = "InlineNS";
		name.Append (fNumNodeSets+1);
		RegisterNodeSet (in, name);

		/* account for no sets or all nodes */
		int index = NodeSetIndex(name);
		if (index > -1)
		{
			/* return name */
			ID.Allocate(1);
			ID[0] = name;
		
			fMessage << " Node Set Name . . . . . . . . . . . . . . . . . = ";
			fMessage << name << '\n';
			fMessage << " Node Set Index. . . . . . . . . . . . . . . . . = ";
			fMessage << index << '\n';
			fMessage << " Node Set Length . . . . . . . . . . . . . . . . = ";
			fMessage << fNodeSetDimensions[index] << '\n';
		}
		else /* empty list */
			ID.Allocate(0);
    }
  else
    {
      int num_sets;
      in >> num_sets;

      ID.Allocate (num_sets);
      for (int i=0; i < num_sets; i++)
	{
	  StringT& name = ID[i];
	  in >> name;
	  int index = NodeSetIndex (name);
	  if (index < 0) {
	  	cout << "\n ModelManagerT::NodeSetList: error retrieving node set " << name << endl;
	  	throw eDatabaseFail;
	  }

	  fMessage << " Node Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Node Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << index << '\n';
	  fMessage << " Node Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fNodeSetDimensions[index] << '\n';
	}
    }
}

void ModelManagerT::SideSetList (ifstreamT& in, ArrayT<StringT>& ID, 
	bool multidatabasesets)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      bool local = true;

      StringT blockID;
      in >> blockID;

      StringT name = "InlineSS";
      name.Append (fNumSideSets+1);
      RegisterSideSet (in, name, local, blockID);

      /* account for no sets */
      int index = SideSetIndex (name);
      if (index > -1)
	{
      ID.Allocate (1);
      ID[0] = name;

	  fMessage << " Side Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Side Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << index << '\n';
	  fMessage << " Side Set Element Group Name . . . . . . . . . . = ";
	  fMessage << blockID << '\n';
	  fMessage << " Side Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fSideSetDimensions[index] << '\n';
	}
	else /* empty list */
		ID.Allocate(0);
    }
  else
    {
      int num_sets;
      if (multidatabasesets)
		in >> num_sets;
      else
		num_sets = 1;

      ID.Allocate (num_sets);
      for (int i=0; i < num_sets; i++)
	{
	  StringT& name = ID[i];
	  in >> name;
	  int index = SideSetIndex (name);
	  if (index < 0) {
	  	cout << "\n ModelManagerT::SideSetList: error retrieving side set " << name << endl;
	  	throw eDatabaseFail;
	  }

	  fMessage << " Side Set Name . . . . . . . . . . . . . . . . . = ";
	  fMessage << name << '\n';
	  fMessage << " Side Set Index. . . . . . . . . . . . . . . . . = ";
	  fMessage << index << '\n';
	  fMessage << " Side Set Element Group Name . . . . . . . . . . = ";
	  fMessage << fInput->SideSetGroupName (name) << '\n';
	  fMessage << " Side Set Length . . . . . . . . . . . . . . . . = ";
	  fMessage << fSideSetDimensions[index] << '\n';
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

//external file needs to opened before ReadCards is called
#if 0
  /* account for text file name instead of data */
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, out, true, "ModelManagerT::ReadCards: could not open file");
#endif

  int count = 0;
  int *pd = data.Pointer();
  StringT ID;
  for (int i=0; i < numc; i++)
    {
      /* read node set name or node number */
      in >> ID;

      /* read nodes in set or create a set from node number */
      if (fFormat == IOBaseT::kTahoe)
	{
	  nodes[i].Allocate (1);
	  nodes[i] = atoi (ID) - 1; // offset
	  count++;
	}
      else
	{
	  nodes[i] = NodeSet (ID);
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
	in >> *pd++;
      in >> value[i];
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

void ModelManagerT::ReadTractionSetData (ifstreamT& in, StringT& element_ID, int& setsize)
{
  if (fFormat == IOBaseT::kTahoe)
    in >> element_ID >> setsize;
  else
    {
      setsize = 1;
      // blockindex set later
    }
}

void ModelManagerT::ReadTractionSideSet (ifstreamT& in, StringT& element_ID, iArray2DT& localsides)
{
  if (fFormat == IOBaseT::kTahoe)
    {
      localsides.Allocate (2,1);
      in >> localsides[0] >> localsides[1];
      // blockindex is already set
    }
  else
    {
		/* read set */
		StringT ss_ID;
		in >> ss_ID; 
		localsides = SideSet(ss_ID);

		/* non-empty set */
		if (localsides.MajorDim() > 0)
		{
			/* try to resolve associated element group ID */
			element_ID = SideSetGroupID(ss_ID);

			/* this shouldn't happen */
			if (!IsSideSetLocal(ss_ID))
			{
				iArray2DT temp = localsides;
				SideSetGlobalToLocal(temp, localsides, element_ID);
			}
		}
		else
			element_ID.Clear();

      fMessage << " Database side set name. . . . . . . . . . . . . = ";
      fMessage << ss_ID << '\n';
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

int ModelManagerT::ElementGroupIndex (const StringT& ID) const
{
  // account for space padding at end of name
  int length1 = ID.Length();
  for (int i=0; i < fNumElementSets; i++)
    {
      int length2 = fElementNames[i].Length();
      int length = (length1 < length2) ? length1 : length2;
      if (strncmp (ID.Pointer(), fElementNames[i].Pointer(), length-1) == 0)
	return i;
    }
  return -1;
}

void ModelManagerT::ElementGroupDimensions (const StringT& ID, int& numelems, int& numelemnodes) const
{
	int index = ElementGroupIndex(ID);
	if (index == -1) {
		cout << "\n ModelManagerT::ElementGroupDimensions: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
	numelems = fElementLengths[index];
	numelemnodes = fElementNodes[index];

//why accept a bad index?
#if 0
  numelems = -1;
  numelemnodes = -1;
  if (index == -1)
    return;
  numelems = fElementLengths[index];
  numelemnodes = fElementNodes[index];
#endif
}

GeometryT::CodeT ModelManagerT::ElementGroupGeometry (const StringT& ID) const
{
	int index = ElementGroupIndex(ID);
	if (index == -1) {
		cout << "\n ModelManagerT::ElementGroupGeometry: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
	return fElementCodes[index];	

//why accept a bad index?
#if 0
  if (index == -1)
    return GeometryT::kNone;
  return fElementCodes[index];
#endif
}

const iArray2DT& ModelManagerT::ElementGroup (const StringT& ID)
{
	ReadConnectivity (ID);
	int index = ElementGroupIndex(ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::ElementGroup: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
	return fElementSets [index];
}

void ModelManagerT::ReadConnectivity (const StringT& ID)
{
	int index = ElementGroupIndex(ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::ReadConnectivity: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
	
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
      fInput->ReadConnectivity (ID, fElementSets[index]);
    }
}

const iArray2DT* ModelManagerT::ElementGroupPointer (const StringT& ID) const
{
	int index = ElementGroupIndex(ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::ElementGroupPointer: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
	return fElementSets.Pointer(index);
}

void ModelManagerT::AllNodeMap (iArrayT& map)
{
	/* no input for kTahoe format */
	if (fFormat == IOBaseT::kTahoe)
	{
		if (map.Length() != fCoordinates.MajorDim()) {
    		cout << "\n ModelManagerT::AllNodeMap: map array is length " << map.Length()
    	         << ", expecting length " << fCoordinates.MajorDim() << endl;
			throw eSizeMismatch;	
		}

		/* default map */
		map.SetValueToPosition();
	}
	else
	{
		if (!fInput) {
    		cout << "\n ModelManagerT::AllNodeMap: input source is not initialized" << endl;
			throw eDatabaseFail;
		}
		if (map.Length() != fInput->NumNodes()) {
    		cout << "\n ModelManagerT::AllNodeMap: map array is length " << map.Length()
    	         << ", expecting length " << fInput->NumNodes() << endl;
			throw eSizeMismatch;	
		}
		fInput->ReadNodeMap (map);
	}
}

void ModelManagerT::AllElementMap (iArrayT& map)
{
	/* no input for kTahoe format */
	if (fFormat == IOBaseT::kTahoe)
	{
		int num_elem = 0;
		for (int i = 0; i < fElementSets.Length(); i++)
			num_elem += fElementSets[i].MajorDim();
		if (map.Length() != num_elem) {
	    	cout << "\n ModelManagerT::AllElementMap: map array is length " << map.Length()
	             << ", expecting length " << num_elem << endl;
			throw eSizeMismatch;	
		}
	
		/* default map */
		map.SetValueToPosition();
	}
	else
	{
		if (!fInput) {
    		cout << "\n ModelManagerT::AllElementMap: input source is not initialized" << endl;
			throw eDatabaseFail;
		}
		if (map.Length() != fInput->NumGlobalElements()) {
	    	cout << "\n ModelManagerT::AllElementMap: map array is length " << map.Length()
	             << ", expecting length " << fInput->NumGlobalElements() << endl;
			throw eSizeMismatch;	
		}
		fInput->ReadAllElementMap (map);
	}
}

void ModelManagerT::ElementMap (const StringT& ID, iArrayT& map)
{
	/* no input for kTahoe format */
	if (fFormat == IOBaseT::kTahoe)
	{
		const iArray2DT& connects = ElementGroup(ID);
		if (map.Length() != connects.MajorDim()) {
    		cout << "\n ModelManagerT::ElementMap: map array is length " << map.Length()
    	         << ", expecting length " << connects.MajorDim() << endl;
			throw eSizeMismatch;		
		}
		
		/* default map */
		map.SetValueToPosition();
	}
	else
	{
		if (!fInput) {
    		cout << "\n ModelManagerT::ElementMap: input source is not initialized" << endl;
			throw eDatabaseFail;
		}
		if (map.Length() != fInput->NumElements(ID)) {
    		cout << "\n ModelManagerT::ElementMap: map array is length " << map.Length()
    	         << ", expecting length " << fInput->NumElements(ID) << endl;
			throw eSizeMismatch;	
		}
		fInput->ReadGlobalElementMap (ID, map);
	}
}

int ModelManagerT::NodeSetIndex (const StringT& ID) const
{
  // account for space padding at end of name
  int length1 = ID.Length();
  for (int i=0; i < fNumNodeSets; i++)
    {
      int length2 = fNodeSetNames[i].Length();
      int length = (length1 < length2) ? length1 : length2;
      if (strncmp (ID.Pointer(), fNodeSetNames[i].Pointer(), length-1) == 0)
	return i;
    }
  return -1;
}

int ModelManagerT::NodeSetLength (const StringT& ID) const
{
	int index = NodeSetIndex(ID);	
	if (index == -1) {
		cout << "\n ModelManagerT::NodeSetLength: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
	
//    return -1; why allow bad name?
  return fNodeSetDimensions [index];
}

const iArrayT& ModelManagerT::NodeSet (const StringT& ID)
{
	int index = NodeSetIndex(ID);	
	if (index == -1) {
		cout << "\n ModelManagerT::NodeSet: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}

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
      fInput->ReadNodeSet (ID, fNodeSets[index]);
    }
  return fNodeSets [index];
}

void ModelManagerT::ManyNodeSets (const ArrayT<StringT>& ID, iArrayT& nodes)
{
  iAutoArrayT temp;
  iArrayT tn;
  for (int i=0; i < ID.Length(); i++)
    {
      tn = NodeSet (ID[i]);
      temp.AppendUnique(tn);
    }

  nodes.Allocate (temp.Length());
  nodes.CopyPart (0, temp, 0, temp.Length());
  nodes.SortAscending ();
}

int ModelManagerT::SideSetIndex (const StringT& ID) const
{
  // account for space padding at end of name
  int length1 = ID.Length();
  for (int i=0; i < fNumSideSets; i++)
    {
      int length2 = fSideSetNames[i].Length();
      int length = (length1 < length2) ? length1 : length2;
      if (strncmp (ID.Pointer(), fSideSetNames[i].Pointer(), length-1) == 0)
	return i;
    }
  return -1;
}

int ModelManagerT::SideSetLength (const StringT& ID) const
{
	int index = SideSetIndex(ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::SideSetLength: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
	return fSideSetDimensions [index];
}

const iArray2DT& ModelManagerT::SideSet (const StringT& ID) const
{
	int index = SideSetIndex(ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::SideSet: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }

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

bool ModelManagerT::IsSideSetLocal (const StringT& ID) const
{
	int index = SideSetIndex(ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::IsSideSetLocal: ID not found: " << ID << endl;
    	throw eOutOfRange;
    }
	return fSideSetIsLocal [index];
}

const StringT& ModelManagerT::SideSetGroupID (const StringT& ss_ID) const
{
	int index = SideSetIndex(ss_ID);
	if (index == -1) {
    	cout << "\n ModelManagerT::SideSetGroupID: ID not found: " << ss_ID << endl;
    	throw eOutOfRange;
    }

	int ss_group_index = fSideSetGroupIndex[index];
	if (ss_group_index < 0 || ss_group_index >= fElementNames.Length()) {
		cout << "\n ModelManagerT::SideSetGroupID: group ID for not defined for set " 
		     << ss_ID << endl;
		throw eOutOfRange;
	} 
	return fElementNames[ss_group_index];
	
	/* need to check if group index < 0, then have global side set and
       need to determine correct group index */
       
	/* group index not defined if side set is empty */
}

void ModelManagerT::SideSetLocalToGlobal (const StringT& element_ID, const iArray2DT& local, iArray2DT& global)
{
	int localelemindex = ElementGroupIndex(element_ID);
	if (localelemindex == -1) {
		cout << "\n ModelManagerT::SideSetLocalToGlobal: element ID not found " << element_ID << endl;
		throw eOutOfRange;
	}

  int offset = 0;
  for (int i=0; i < localelemindex; i++)
    offset += fElementLengths[i];

  global = local;
  int *pelem = global.Pointer();
  for (int j=0; j < global.MajorDim(); j++, pelem += 2)
    *pelem += offset;
}

void ModelManagerT::SideSetGlobalToLocal(const iArray2DT& global, iArray2DT& local, 
	StringT& element_ID)
{
#pragma unused(element_ID)
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

bool ModelManagerT::RegisterVariElements (const StringT& ID, nVariArray2DT<int>& conn, 
					  GeometryT::CodeT code, int numelemnodes,
					  int headroom)
{
  if (!CheckID (fElementNames, ID, "Element Group")) return false;
  
  fElementNames.Append (ID);
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
void ModelManagerT::UpdateConnectivity (const StringT& ID, const iArray2DT& connects)
{
	int index = ElementGroupIndex(ID);
	if (index == -1) {
		cout << "\n ModelManagerT::UpdateConnectivity: element ID not found " << ID << endl;
		throw eOutOfRange;
	}
	fElementSets[index] = connects;
	fElementLengths[index] = fElementSets[index].MajorDim();
	fElementNodes[index] = fElementSets[index].MinorDim();
}

void ModelManagerT::AddElement (const StringT& ID, const iArray2DT& connects, 
	iArrayT& new_elem_tags, int& numelems)
{
	int index = ElementGroupIndex(ID);
	if (index == -1) {
		cout << "\n ModelManagerT::AddElement: element ID not found " << ID << endl;
		throw eOutOfRange;
	}

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

/*************************************************************************
* Private
*************************************************************************/

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
#ifdef __ACCESS__
    case IOBaseT::kExodusII:
      fInput = new ExodusInputT (fMessage);
      break;
#endif
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

  		/* gather side set info */
		for (int i=0; i < fNumSideSets; i++)
		{
			fSideSetDimensions[i] = fInput->NumSidesInSet (fSideSetNames[i]);
			if (fSideSetIsLocal[i] && 
			    fSideSetDimensions[i] > 0) /* don't try this with an empty set */
			{
				const StringT& name = fInput->SideSetGroupName(fSideSetNames[i]);
				fSideSetGroupIndex[i] = ElementGroupIndex(name);
			}
		}
    }
  return true;
}

bool ModelManagerT::CheckID (const ArrayT<StringT>& list, const StringT& ID, const char *type) const
{
  // account for space padding at end of name
  int l1 = ID.Length();

  for (int i=0; i < list.Length(); i++)
    {
      int l2 = list[i].Length();
      int l = (l1 < l2) ? l1 : l2;
      if (strncmp (list[i].Pointer(), ID.Pointer(), l-1) == 0)
	{
	  fMessage << "\nModelManagerT::CheckID\n";
	  fMessage << "   " << type << " already has a registered set called " << ID << "\n\n";
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
