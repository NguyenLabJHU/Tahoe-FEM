/* $Id: ModelManagerT.cpp,v 1.24 2002-03-26 17:24:08 paklein Exp $ */
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

/* array behavior */
const bool ArrayT<ModelManagerT::SideSetScopeT>::fByteCopy = true;

ModelManagerT::ModelManagerT (ostream& message):
	fMessage(message),
	fCoordinateDimensions (2),
	fInput(NULL)
{
  fCoordinateDimensions = -1;
}

ModelManagerT::~ModelManagerT(void) { Clear(); }

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

bool ModelManagerT::RegisterElementGroup (ifstreamT& in, const StringT& ID, GeometryT::CodeT code)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterElementGroup(ifstreamT): count not open file");

  int length, numnodes;
  in2 >> length >> numnodes;
  iArray2DT temp (length, numnodes);
  temp.ReadNumbered (in2);
  temp += -1;
  return RegisterElementGroup (ID, temp, code, true);
}

bool ModelManagerT::RegisterElementGroup (const StringT& ID, iArray2DT& conn, 
	GeometryT::CodeT code, bool keep)
{
	if (!CheckID (fElementNames, ID, "Element Group")) return false;
  
  	/* element group parameters */
	fElementNames.Append (ID);
	fElementLengths.Append (conn.MajorDim());
	fElementNodes.Append (conn.MinorDim());
	fElementCodes.Append (code);

	iArray2DT* set = NULL;
	if (!keep || !conn.IsAllocated()) /* make copy */
		set = new iArray2DT(conn);
	else /* take memory */
	{
		set = new iArray2DT;
		set->Swap(conn);
		conn.Alias(*set);
	}
	fElementSets.Append(set);

	return true;
}

bool ModelManagerT::RegisterNodeSet (const StringT& ID, iArrayT& set, bool keep)
{
	if (!CheckID (fNodeSetNames, ID, "Node Set")) return false;
  
	fNodeSetNames.Append (ID);
	fNodeSetDimensions.Append (set.Length());

	iArrayT* new_set = NULL;
	if (!keep || !set.IsAllocated()) /* make copy */
		new_set = new iArrayT(set);
	else /* take memory */
	{
		new_set = new iArrayT;
		new_set->Swap(set);
		set.Alias(*new_set);
	}
  	fNodeSets.Append(new_set);
  	return true;
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
      return RegisterNodeSet(ID, n, true);
    }
  else
    return false;
}

void ModelManagerT::ReadInlineCoordinates (ifstreamT& in)
{
  if (fFormat == IOBaseT::kTahoe)
    RegisterNodes (in);
}

bool ModelManagerT::RegisterNodes (ifstreamT& in)
{
  ifstreamT tmp;
  ifstreamT& in2 = OpenExternal (in, tmp, fMessage, true, "ModelManagerT::RegisterNodes(ifstreamT): count not open file");

  in2 >> fCoordinateDimensions[0] >> fCoordinateDimensions[1];
  fCoordinates.Dimension(fCoordinateDimensions[0], fCoordinateDimensions[1]);
  fCoordinates.ReadNumbered (in2);
  return true;
}

bool ModelManagerT::RegisterNodes(dArray2DT& coords, bool keep)
{
	if (!keep || !coords.IsAllocated())
		fCoordinates = coords;
	else
	{
		fCoordinates.Swap(coords);
		coords.Alias(fCoordinates);
	}
	
	fCoordinateDimensions[0] = fCoordinates.MajorDim();
	fCoordinateDimensions[1] = fCoordinates.MinorDim();

	return true;
}

bool ModelManagerT::RegisterSideSet (const StringT& ss_ID, iArray2DT& set, 
	SideSetScopeT scope, const StringT& element_ID, bool keep)
{
	if (!CheckID (fSideSetNames, ss_ID, "Side Set")) return false;

	/* side set parameters */
 	fSideSetNames.Append(ss_ID);
	fSideSetDimensions.Append(set.MajorDim());
	fSideSetScope.Append(scope);
 	if (scope == kLocal)
	{
		int index = ElementGroupIndex(element_ID);
		if (index == kNotFound) {
			cout << "\n ModelManagerT::RegisterSideSet: element ID not found: " << element_ID << endl;
			throw eOutOfRange;
		}
    	fSideSetGroupIndex.Append (index);
  	}
  	else
		fSideSetGroupIndex.Append (kNotFound);

	iArray2DT* new_set = NULL;
	if (!keep || !set.IsAllocated()) /* make copy */
		new_set = new iArray2DT(set);
	else /* take memory */
	{
		new_set = new iArray2DT;
		new_set->Swap(set);
		set.Alias(*new_set);
	}
	fSideSets.Append(new_set);
	return true;
}

bool ModelManagerT::RegisterSideSet (ifstreamT& in, const StringT& ss_ID, SideSetScopeT scope, 
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
      return RegisterSideSet (ss_ID, s, scope, element_ID, true);
    }
  else
    return false;
}

void ModelManagerT::ElementBlockList (ifstreamT& in, ArrayT<StringT>& ID, iArrayT& matnums)
{
  /* number of blocks in element data */
  int num_blocks = 0;
  in >> num_blocks;
  fMessage << " Number of connectivity data blocks. . . . . . . = " << num_blocks << '\n';
  if (num_blocks < 1) throw eBadInputValue;

  ID.Dimension (num_blocks);
  matnums.Dimension (num_blocks);
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
	  name.Append (NumElementGroups() + 1);
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
		name.Append (NumNodeSets() + 1);
		RegisterNodeSet (in, name);

		/* account for no sets or all nodes */
		int index = NodeSetIndex(name);
		if (index > kNotFound)
		{
			/* return name */
			ID.Dimension(1);
			ID[0] = name;
		
			fMessage << " Node Set Name . . . . . . . . . . . . . . . . . = ";
			fMessage << name << '\n';
			fMessage << " Node Set Index. . . . . . . . . . . . . . . . . = ";
			fMessage << index << '\n';
			fMessage << " Node Set Length . . . . . . . . . . . . . . . . = ";
			fMessage << fNodeSetDimensions[index] << '\n';
		}
		else /* empty list */
			ID.Dimension(0);
    }
  else
    {
      int num_sets;
      in >> num_sets;

      ID.Dimension (num_sets);
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
      StringT blockID;
      in >> blockID;

      StringT name = "InlineSS";
      name.Append (NumSideSets() + 1);
      RegisterSideSet (in, name, kLocal, blockID);

      /* account for no sets */
      int index = SideSetIndex (name);
      if (index > kNotFound)
	{
      ID.Dimension (1);
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
		ID.Dimension(0);
    }
  else
    {
      int num_sets;
      if (multidatabasesets)
		in >> num_sets;
      else
		num_sets = 1;

      ID.Dimension (num_sets);
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
	  fMessage << Input().SideSetGroupName (name) << '\n';
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
	  nodes[i].Dimension (1);
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
      localsides.Dimension (2,1);
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
	ReadCoordinates ();
	return fCoordinates;
}

void ModelManagerT::ReadCoordinates(void)
{
	/* not yet loaded */
	if (fCoordinates.MajorDim() != fCoordinateDimensions[0] ||
	    fCoordinates.MajorDim() != fCoordinateDimensions[1])
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
	
		fCoordinates.Dimension(fCoordinateDimensions[0], fCoordinateDimensions[1]); 
		Input("ReadCoordinates").ReadCoordinates (fCoordinates);
	}
}

/* used to reduce 3D database coordinates (Patran, Abaqus, etc.) */
bool ModelManagerT::AreElements2D (void) const
{
  if (fCoordinateDimensions[1] < 3) return true;

  /* look over registered element sets */
  for (int i=0; i < NumElementGroups(); i++)
    if (fElementCodes[i] == GeometryT::kPoint ||
	fElementCodes[i] == GeometryT::kTriangle ||
	fElementCodes[i] == GeometryT::kQuadrilateral )
      return true;

  return false;
}

/* return an unused element ID */
StringT ModelManagerT::FreeElementID(const StringT& prefix) const
{
	int tail = 0;
	bool free = false;
	StringT ID;
	while (!free)
	{
		ID.Clear();
		ID.Append(prefix, tail);
		free = true;
		for (int i = 0; free && i < fElementNames.Length(); i++)
			if (ID == fElementNames[i])
				free = false;
		tail++;
	}
	return ID;
}

int ModelManagerT::ElementGroupIndex (const StringT& ID) const
{
	/* scan element names */
	int num_sets = NumElementGroups();
  	for (int i=0; i < num_sets; i++)
		if (ID_Match(ID, fElementNames[i]))
			return i;
	
	/* not found */
	return kNotFound;
}

void ModelManagerT::ElementGroupDimensions (const StringT& ID, int& numelems, int& numelemnodes) const
{
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
		cout << "\n ModelManagerT::ElementGroupDimensions: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
	numelems = fElementLengths[index];
	numelemnodes = fElementNodes[index];

//why accept a bad index?
#if 0
  numelems = -1;
  numelemnodes = -1;
  if (index == kNotFound)
    return;
  numelems = fElementLengths[index];
  numelemnodes = fElementNodes[index];
#endif
}

GeometryT::CodeT ModelManagerT::ElementGroupGeometry (const StringT& ID) const
{
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
		cout << "\n ModelManagerT::ElementGroupGeometry: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
	return fElementCodes[index];	

//why accept a bad index?
#if 0
  if (index == kNotFound)
    return GeometryT::kNone;
  return fElementCodes[index];
#endif
}

const iArray2DT& ModelManagerT::ElementGroup (const StringT& ID)
{
	ReadConnectivity (ID);
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
    	cout << "\n ModelManagerT::ElementGroup: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
    
    iArray2DT* set = fElementSets[index];
    if (!set) {
    	cout << "\n ModelManagerT::ElementGroup: set " << ID 
    	     << " not initialized" << endl;
    }
    
	return *set;
}

void ModelManagerT::ReadConnectivity (const StringT& ID)
{
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
    	cout << "\n ModelManagerT::ReadConnectivity: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
    if (!fElementSets[index]) {
    	cout << "\n ModelManagerT::ReadConnectivity: set " << ID 
    	     << " is not initialized" << endl;    
    	throw eGeneralFail;
    }
	
	/* data not yet loaded */
	if (fElementSets[index]->MajorDim() != fElementLengths[index] ||
	    fElementSets[index]->MinorDim() != fElementNodes[index])
	{
		if (fFormat == IOBaseT::kTahoe)
		{
			cout << "\n ModelManagerT::ReadConnectivity, elems not registered yet" << endl;
			throw eGeneralFail;
		}
		
		/* allocate space */
		fElementSets[index]->Dimension(fElementLengths[index], fElementNodes[index]);

		/* do read */
		try { Input("ReadConnectivity").ReadConnectivity (ID, *fElementSets[index]); }
		catch (int exception) {
			cout << "\n ModelManagerT::ReadConnectivity: exception reading ID " << ID 
			     << ": " << exception << endl;
			throw eDatabaseFail;
		}
    }
}

const iArray2DT* ModelManagerT::ElementGroupPointer (const StringT& ID) const
{
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
    	cout << "\n ModelManagerT::ElementGroupPointer: ID not found: " << ID << endl;
    	throw eOutOfRange;
    }
	return fElementSets[index];
}

void ModelManagerT::AllNodeIDs (iArrayT& ids)
{
	/* no input for kTahoe format */
	if (fFormat == IOBaseT::kTahoe)
	{
		if (ids.Length() != fCoordinates.MajorDim()) {
    		cout << "\n ModelManagerT::AllNodeIDs: ids array is length " << ids.Length()
    	         << ", expecting length " << fCoordinates.MajorDim() << endl;
			throw eSizeMismatch;	
		}

		/* default ids */
		ids.SetValueToPosition();
	}
	else
	{
		InputBaseT& input = Input("AllNodeIDs");

		/* dimension (check) */
		if (ids.Length() == 0)
			ids.Dimension(input.NumNodes());		
		else if (ids.Length() != input.NumNodes()) {
    		cout << "\n ModelManagerT::AllNodeIDs: ids array is length " << ids.Length()
    	         << ", expecting length " << input.NumNodes() << endl;
			throw eSizeMismatch;	
		}
		
		/* read */
		input.ReadNodeID(ids);
	}
}

void ModelManagerT::AllElementIDs (iArrayT& ids)
{
	/* no input for kTahoe format */
	if (fFormat == IOBaseT::kTahoe)
	{
		int num_elem = 0;
		for (int i = 0; i < fElementSets.Length(); i++)
			num_elem += fElementSets[i]->MajorDim();
		if (ids.Length() != num_elem) {
	    	cout << "\n ModelManagerT::AllElementIDs: ids array is length " << ids.Length()
	             << ", expecting length " << num_elem << endl;
			throw eSizeMismatch;	
		}
	
		/* default ids */
		ids.SetValueToPosition();
	}
	else
	{
		InputBaseT& input = Input("AllElementIDs");
		if (ids.Length() != input.NumGlobalElements()) {
	    	cout << "\n ModelManagerT::AllElementIDs: ids array is length " << ids.Length()
	             << ", expecting length " << input.NumGlobalElements() << endl;
			throw eSizeMismatch;	
		}
		input.ReadAllElementMap (ids);
	}
}

void ModelManagerT::ElementIDs (const StringT& ID, iArrayT& ids)
{
	/* no input for kTahoe format */
	if (fFormat == IOBaseT::kTahoe)
	{
		const iArray2DT& connects = ElementGroup(ID);
		if (ids.Length() != connects.MajorDim()) {
    		cout << "\n ModelManagerT::ElementIDs: ids array is length " << ids.Length()
    	         << ", expecting length " << connects.MajorDim() << " for ID " << ID << endl;
			throw eSizeMismatch;		
		}
		
		/* default ids */
		ids.SetValueToPosition();
	}
	else
	{
		InputBaseT& input = Input("ElementIDs");
		if (ids.Length() != input.NumElements(ID)) {
    		cout << "\n ModelManagerT::ElementIDs: ids array is length " << ids.Length()
    	         << ", expecting length " << input.NumElements(ID) << " for ID " 
    	         << ID << endl;
			throw eSizeMismatch;	
		}
		input.ReadGlobalElementMap (ID, ids);
	}
}

void ModelManagerT::AllNodeMap (iArrayT& map)
{
	if (map.Length() != fCoordinateDimensions[0])
	{
		cout << "\n ModelManagerT::NodeMap: map array is length " << map.Length()
             << ", expecting length " << fCoordinates.MajorDim() << endl;
		throw eSizeMismatch;	
	}
	map.SetValueToPosition();
}

void ModelManagerT::ElementMap (const StringT& ID, iArrayT& map)
{
  /* no input for kTahoe format */
  if (fFormat == IOBaseT::kTahoe)
    {
      const iArray2DT& connects = ElementGroup (ID);
      if (map.Length() != connects.MajorDim()) {
	cout << "\n ModelManagerT::ElementMap: map array is length " << map.Length()
	     << ", expecting length " << connects.MajorDim() << " for ID " << ID << endl;
	throw eSizeMismatch;	
      }
      /* default map */
      map.SetValueToPosition();
    }
	else
	{
		InputBaseT& input = Input ("Element Set");
		input.ReadGlobalElementSet (ID, map);
  	}
}

int ModelManagerT::NodeSetIndex (const StringT& ID) const
{
	/* scan node set names */
	int num_sets = NumNodeSets(); 
  	for (int i = 0; i < num_sets; i++)
  		if (ID_Match(ID, fNodeSetNames[i]))
  			return i;
  			
  	/* not found */
  	return kNotFound;
}

int ModelManagerT::NodeSetLength (const StringT& ID) const
{
	int index = NodeSetIndex(ID);	
	if (index == kNotFound) {
		cout << "\n ModelManagerT::NodeSetLength: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
	
//    return -1; why allow bad name?
  return fNodeSetDimensions [index];
}

const iArrayT& ModelManagerT::NodeSet (const StringT& ID)
{
	int index = NodeSetIndex(ID);	
	if (index == kNotFound) {
		cout << "\n ModelManagerT::NodeSet: ID not found: " << ID << endl;
		throw eDatabaseFail;
	}
    if (!fNodeSets[index]) {
    	cout << "\n ModelManagerT::NodeSet: set " << ID 
    	     << " is not initialized" << endl;    
    	throw eGeneralFail;
    }

	/* not yet loaded */
	if (fNodeSets[index]->Length() != fNodeSetDimensions[index])
	{
		if (fFormat == IOBaseT::kTahoe)
		{
			cout << "\n ModelManagerT::NodeSet, set not registered yet" << endl;
			throw eGeneralFail;
		}
		fNodeSets[index]->Dimension (fNodeSetDimensions[index]);
		try { Input("NodeSet").ReadNodeSet(ID, *fNodeSets[index]); }
		catch (int exception) {
			cout << "\n ModelManagerT::NodeSet: exception reading ID " << ID 
			     << ": " << exception << endl;
			throw eDatabaseFail;
		}
	}
	return *fNodeSets[index];
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

  nodes.Dimension (temp.Length());
  nodes.CopyPart (0, temp, 0, temp.Length());
  nodes.SortAscending ();
}

int ModelManagerT::SideSetIndex (const StringT& ID) const
{
	/* scan side set names */
	int num_sets = NumSideSets();
	for (int i = 0; i < num_sets; i++)
		if (ID_Match(ID, fSideSetNames[i]))
			return i;

	/* not found */
	return kNotFound;
}

int ModelManagerT::SideSetLength (const StringT& ID) const
{
	int index = SideSetIndex(ID);
	if (index == kNotFound) {
    	cout << "\n ModelManagerT::SideSetLength: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
	return fSideSetDimensions [index];
}

const iArray2DT& ModelManagerT::SideSet (const StringT& ID) const
{
	int index = SideSetIndex(ID);
	if (index == kNotFound) {
    	cout << "\n ModelManagerT::SideSet: ID not found: " << ID<< endl;
    	throw eOutOfRange;
    }
    if (!fSideSets[index]) {
    	cout << "\n ModelManagerT::SideSet: set " << ID 
    	     << " is not initialized" << endl;    
    	throw eGeneralFail;
    }

	if (fSideSets[index]->MajorDim() != fSideSetDimensions[index])
    {
		if (fFormat == IOBaseT::kTahoe)
		{
			cout << "\n ModelManagerT::SideSet, set not registered yet" << endl;
			throw eGeneralFail;
		}
		InputBaseT& input = Input("SideSet");
		fSideSets[index]->Dimension (fSideSetDimensions[index], 2);
		try {
			if (fSideSetScope[index] == kLocal)
				input.ReadSideSetLocal (fSideSetNames[index], *fSideSets[index]);
			else
				input.ReadSideSetGlobal (fSideSetNames[index], *fSideSets[index]);
		}
		catch (int exception) {
			cout << "\n ModelManagerT::SideSet: exception reading ID " << ID 
			     << ": " << exception << endl;
			throw eDatabaseFail;
		}
	}
	return *fSideSets[index];
}

bool ModelManagerT::IsSideSetLocal (const StringT& ID) const
{
	int index = SideSetIndex(ID);
	if (index == kNotFound) {
    	cout << "\n ModelManagerT::IsSideSetLocal: ID not found: " << ID << endl;
    	throw eOutOfRange;
    }
	return fSideSetScope[index] == kLocal;
}

const StringT& ModelManagerT::SideSetGroupID (const StringT& ss_ID) const
{
	int index = SideSetIndex(ss_ID);
	if (index == kNotFound) {
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
	if (localelemindex == kNotFound) {
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
  RegisterNodes (temp, true);
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
  iArray2DT* set = new iArray2DT;
  fElementSets.Append(set);

  conn.SetWard (headroom, *fElementSets[index], numelemnodes);
  return true;
}

/* call this function after the connectivity has been changed by outside classes */
void ModelManagerT::UpdateConnectivity (const StringT& ID, iArray2DT& connects, bool keep)
{
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
		cout << "\n ModelManagerT::UpdateConnectivity: element ID not found " << ID << endl;
		throw eOutOfRange;
	}
	if (!fElementSets[index]) throw eGeneralFail;
	
	if (!keep || !connects.IsAllocated())
		*fElementSets[index] = connects;
	else
	{
		fElementSets[index]->Swap(connects);
		connects.Alias(*fElementSets[index]);
	}
	fElementLengths[index] = fElementSets[index]->MajorDim();
	fElementNodes[index] = fElementSets[index]->MinorDim();
}

void ModelManagerT::AddElement (const StringT& ID, const iArray2DT& connects, 
	iArrayT& new_elem_tags, int& numelems)
{
	int index = ElementGroupIndex(ID);
	if (index == kNotFound) {
		cout << "\n ModelManagerT::AddElement: element ID not found " << ID << endl;
		throw eOutOfRange;
	}
	if (!fElementSets[index]) throw eGeneralFail;

  if (connects.MajorDim() != new_elem_tags.Length() ||
      connects.MinorDim() != fElementSets[index]->MinorDim() ) throw eSizeMismatch;

  /* set new elem tags to the old last elem number + 1 */
  new_elem_tags.SetValueToPosition ();
  new_elem_tags += fElementSets[index]->MajorDim();

  /* reset the number of elements */
  int newelems = connects.MajorDim();
  fElementLengths[index] += newelems;
  numelems = fElementLengths[index];

  /* reallocate */
  fElementSets[index]->Resize (newelems);

  /* copy in */
  int *pc = connects.Pointer();
  for (int i=0; i < newelems; i++, pc += connects.MinorDim())
    fElementSets[index]->SetRow (new_elem_tags[i], pc);
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

/** return true of the ID's match */
bool ModelManagerT::ID_Match(const StringT& a, const StringT& b) const
{
	int la = ID_Length(a);
	int lb = ID_Length(b);
	if (la != lb)
		return false;
	else
		return strncmp(a, b, la) == 0;
}

/** return the length of the ID string not including any trailing white-space */
int ModelManagerT::ID_Length(const StringT& ID) const
{
	const char* str = ID.Pointer();
	int count = 0;
	int length = strlen(str);
	while (count < length && !isspace(*str++)) 
		count++;
	return count;
}

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
	
  int num_elem_sets = fInput->NumElementGroups ();
  fElementLengths.Dimension (num_elem_sets);
  fElementNodes.Dimension (num_elem_sets);
  fElementNames.Dimension (num_elem_sets);
  fElementCodes.Dimension (num_elem_sets);
  fElementSets.Dimension (num_elem_sets);
  fElementSets = NULL;
	if (num_elem_sets > 0)
	{
		fInput->ElementGroupNames (fElementNames);
		for (int e=0; e < num_elem_sets; e++)
		{
			fElementLengths[e] = fInput->NumElements (fElementNames[e]);
			fElementNodes[e] = fInput->NumElementNodes (fElementNames[e]);
			fInput->ReadGeometryCode (fElementNames[e], fElementCodes[e]);
			
			/* create empty set */
			fElementSets[e] = new iArray2DT;
		}
	}
	return true;
}

bool ModelManagerT::ScanNodeSets (void)
{
	/* check if input has been initialized */
	if (!fInput) return false;

  int num_node_sets = fInput->NumNodeSets();
  fNodeSetNames.Dimension (num_node_sets);
  fNodeSetDimensions.Dimension (num_node_sets);
  fNodeSets.Dimension (num_node_sets);
  fNodeSets = NULL;
	if (num_node_sets > 0)
    {
		fInput->NodeSetNames(fNodeSetNames);
		for (int i=0; i < num_node_sets; i++)
		{
			fNodeSetDimensions[i] = fInput->NumNodesInSet (fNodeSetNames[i]);
			
			/* create empty set */
			fNodeSets[i] = new iArrayT;
		}
    }
 	 return true;
}

bool ModelManagerT::ScanSideSets (void)
{
	/* check if input has been initialized */
	if (!fInput) return false;

  int num_side_sets = fInput->NumSideSets();
  fSideSetNames.Dimension (num_side_sets);
  fSideSetDimensions.Dimension (num_side_sets);
  fSideSetScope.Dimension (num_side_sets);
  fSideSetGroupIndex.Dimension (num_side_sets);
  fSideSets.Dimension (num_side_sets);
  fSideSets = NULL;	
	if (num_side_sets > 0)
	{
		fInput->SideSetNames (fSideSetNames);
		if (fInput->AreSideSetsLocal())
			fSideSetScope = kLocal;
		else
			fSideSetScope = kGlobal;
		fSideSetGroupIndex = kNotFound;

  		/* gather side set info */
		for (int i=0; i < num_side_sets; i++)
		{
			fSideSetDimensions[i] = fInput->NumSidesInSet (fSideSetNames[i]);
			if (fSideSetScope[i] == kLocal && 
			    fSideSetDimensions[i] > 0) /* don't try this with an empty set */
			{
				const StringT& name = fInput->SideSetGroupName(fSideSetNames[i]);
				fSideSetGroupIndex[i] = ElementGroupIndex(name);
			}
			
			/* create empty set */
			fSideSets[i] = new iArray2DT(0,2);
		}
    }
  return true;
}

bool ModelManagerT::CheckID (const ArrayT<StringT>& list, const StringT& ID, const char *type) const
{
	for (int i=0; i < list.Length(); i++)
		if (ID_Match(ID, list[i]))
		{
	  		fMessage << "\n ModelManagerT::CheckID\n";
			fMessage << "   " << type << " already has a registered set called " << ID << "\n\n";
			fMessage << "  Sets: \n";
			for (int j=0; j < list.Length(); j++)
				fMessage << "       " << list[i] << "\n";
			fMessage << "\n";
			return false;
		}

	/* OK */
	return true;
}

/* clear database parameters */
void ModelManagerT::Clear(void)
{
	/* protected attributes */
	fFormat = IOBaseT::kNone;
	delete fInput;
	fInput = NULL;
	fInputName.Clear();
	fCoordinateDimensions = -1;

	fElementLengths.Dimension(0);
	fElementNodes.Dimension(0);
	fNodeSetDimensions.Dimension(0);
	fSideSetDimensions.Dimension(0);
	fSideSetScope.Dimension(0);
	fSideSetGroupIndex.Dimension(0);

	fElementNames.Dimension(0);
	fNodeSetNames.Dimension(0);
	fSideSetNames.Dimension(0);
	fElementCodes.Dimension(0);
	
	for (int i = 0; i < fElementSets.Length(); i++)
		delete fElementSets[i];
	fElementSets.Dimension(0);

	for (int i = 0; i < fNodeSets.Length(); i++)
		delete fNodeSets[i];
	fNodeSets.Dimension(0);

	for (int i = 0; i < fSideSets.Length(); i++)
		delete fSideSets[i];
	fSideSets.Dimension(0);
}
