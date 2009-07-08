/* $Id: TranslateIOManager.cpp,v 1.43 2009-07-08 17:56:27 regueiro Exp $  */
#include "TranslateIOManager.h"

#include "ExceptionT.h"
#include "IOBaseT.h"
#include "OutputSetT.h"
#include "fstreamT.h"

using namespace Tahoe;

TranslateIOManager::TranslateIOManager (ostream& out, istream& in, bool write) :
  fMessage (out),
  fIn (in),
  fWrite (write),
  fModel (out),
  fOutput (NULL),
  fNumNV (0),
  fNumEV (0),
  fNumQV (0),
  fNumTS (0),
  fCoords (0),
  fEcho (false)
{
}
TranslateIOManager::~TranslateIOManager(void)
{
  delete fOutput;
}

void TranslateIOManager::SetEcho (int s, const StringT& file)
{
	fEchoOut.open(file);
	if (!fEchoOut.is_open()) 
		ExceptionT::GeneralFail("TranslateIOManager::SetEcho","Cannot open file %s", file.Pointer());
  	fEcho = true;
  	if (fEcho)
  		fEchoOut << "%\n" << s << endl;
}

void TranslateIOManager::Translate (const StringT& program, const StringT& version, const StringT& title)
{
  SetInput ();
  SetOutput (program, version, title);

  InitializeVariables ();
  if (fNumNV < 1 && fNumEV < 1) 
    {
      fMessage << "\n No variables found, writing geometry only.";
      WriteGeometry ();
      fOutput->WriteGeometry ();
      return;
    }
  
  WriteGeometry ();

  InitializeTime ();
  if (fNumTS < 1)
    {
      fMessage << "\n**** No time steps found, writing geometry only.\n";
      fOutput->WriteGeometry ();
      return;      
    }
  
  TranslateVariables ();
}


/**************** PROTECTED **********************/

void TranslateIOManager::SetInput (void)
{
  IOBaseT::FileTypeT format;
  IOBaseT temp (cout);
  StringT database;
  if (fWrite)
    {
      temp.InputFormats (cout);
      cout << "\n Enter the Model Format Type: ";
    }
  fIn >> format;
  if (format != IOBaseT::kTahoe)
    {
      if (fWrite)
	cout << "\n Enter the Model File Name: ";
      fIn >> database;
      database.ToNativePathName();
    }
  else
    database = "\0";
  if (fModel.Initialize (format, database, true))
    cout << "\n Input Format: " << format << " File: " << database << endl;
  else
    ExceptionT::DatabaseFail ("TranslateIOManager::SetInput", "Unable to initialize model file %s", database.Pointer());

  // echo
  if (fEcho) fEchoOut << format << " " << database << endl;
}

void TranslateIOManager::SetOutput (const StringT& program_name, const StringT& version, const StringT& title)
{
	if (fWrite)
	{
		cout << "\n\n";
		IOBaseT::OutputFormats (cout);
		cout << "\n Enter the Output Format: ";
	}
	IOBaseT::FileTypeT file_type;
	fIn >> file_type;

  if (fWrite)
    cout << "\n Enter the root of the output files: ";
  fIn >> fOutputName;
  cout << "\n Output format: " << file_type << "\n File: " << fOutputName << endl; 
  
  // echo
  if (fEcho) fEchoOut << file_type << " " << fOutputName << endl;

  fOutputName.ToNativePathName();
  fOutputName.Append(".ext"); //trimmed off by fOutput
  
	/* construct output formatter */
	fOutput = IOBaseT::NewOutput(program_name, version, title, fOutputName, file_type, cout);
	if (!fOutput) ExceptionT::OutOfMemory ("TranslateIOManager::SetOutput");
}

void TranslateIOManager::InitializeVariables (void)
{
  fNumNV = fModel.NumNodeVariables ();
  fNumEV = fModel.NumElementVariables ();
  fNumQV = fModel.NumQuadratureVariables ();

  cout << "\n" << setw (10) << fNumNV << " Node Variables\n";
  cout << setw (10) << fNumEV << " Element Variables\n";
  cout << setw (10) << fNumQV << " Quadrature Variables\n";

  fNodeLabels.Dimension (fNumNV);
  fKeepNodeLabels.Dimension (fNumNV);
  fElementLabels.Dimension (fNumEV);
  fKeepElementLabels.Dimension (fNumEV);

  if (fNumNV > 0) {
  	fModel.NodeLabels (fNodeLabels);
  	ReNameLabels("node", fNodeLabels);
  	FilterLabels("node", fKeepNodeLabels, fNodeLabels);
  	fNumNV = fNodeLabels.Length();
  }
  
  if (fNumEV > 0) {
  	fModel.ElementLabels (fElementLabels);
  	ReNameLabels("element", fElementLabels);
  	FilterLabels("element", fKeepElementLabels, fElementLabels);
  	fNumEV = fElementLabels.Length();
  }
}

void TranslateIOManager::InitializeNodeVariables (void)
{
  fNumNV = fModel.NumNodeVariables ();
  cout << "\n" << setw (10) << fNumNV << " Node Variables\n\n";
  fNodeLabels.Dimension (fNumNV);
  if (fNumNV > 0) fModel.NodeLabels (fNodeLabels);

  // query user as to which variables to translate
  VariableQuery (fNodeLabels, fNVUsed);

  StringT answer;
  if (fWrite)
    cout << "\n Do you wish to translate coordinate values (y/n) ? ";
  fIn >> answer;
  if (fEcho) fEchoOut << answer[0] << endl;
  
  if (answer[0] == 'y' || answer[0] == 'Y')
    {
      fCoords = fModel.NumDimensions();
      int numnodes = fModel.NumNodes();
      cout << "\n Adding coordinates to variable values.\n";
    }
  else
    fCoords = 0;
}

void TranslateIOManager::InitializeElementVariables(void)
{
	fNumEV = fModel.NumElementVariables ();
	cout << "\n" << setw (10) << fNumEV << " Element Variables\n\n";
	fElementLabels.Dimension(fNumEV);
	if (fNumEV > 0) fModel.ElementLabels(fElementLabels);

	// query user as to which variables to translate
	VariableQuery(fElementLabels, fEVUsed);
}

void TranslateIOManager::InitializeQuadVariables (void)
{
  fNumQV = fModel.NumQuadratureVariables ();
  cout << "\n" << setw (10) << fNumQV << " Quadrature Variables\n\n";
  fQuadratureLabels.Dimension (fNumQV);
  if (fNumQV > 0) fModel.QuadratureLabels (fQuadratureLabels);

  // query user as to which variables to translate
  if (fNumQV < 1)
    ExceptionT::GeneralFail ("TranslateIOManager::InitializeQuadVariables","No quadrature variables found");
  //cout << fNumQV << " " << fQuadratureLabels[0] <<endl;
  VariableQuery (fQuadratureLabels, fQVUsed);
}

void TranslateIOManager::InitializeElements (int& group, StringT& groupname)
{
  int num = fModel.NumElementGroups ();
  const ArrayT<StringT>& elemsetnames = fModel.ElementGroupIDs();
  if (fWrite)
    {
      cout << "\n";
      for (int h=0; h < num; h++)
	cout << "    " << h+1 << ". " << elemsetnames[h] << "\n";
      cout << "\n You must have one type of element within the group you select.\n";
      cout << " Enter the number of the element group: ";
    }
  fIn >> group;
  if (fEcho) fEchoOut << group << endl;
  if (group < 1 || group > elemsetnames.Length()) 
    ExceptionT::OutOfRange ("TranslateIOManager::InitializeElements","Element Group Index %i is invalid", group);
  else
    cout << "\n Translating element group: " << group << " " << elemsetnames[group-1] << endl;
  group--;
  groupname = elemsetnames[group];
}

void TranslateIOManager::InitializeNodePoints (iArrayT& nodes, iArrayT& index)
{
  int selection;
  if (fWrite)
    {
      cout << "\n One file will be written per node point.\n";
      cout << "1. List of nodes\n";
      cout << "2. Node Set\n";
      cout << "3. Every nth node\n";
      cout << "4. Start node, end node, and increment\n";
      cout << "\n How do you want to define your list of nodes: ";
    }
  fIn >> selection;
  if (fEcho) fEchoOut << selection << endl;

  int numnodes = fModel.NumNodes();
  int numdims = fModel.NumDimensions();
  iArrayT nodeIDs (numnodes);
  fNodeMap.Dimension (numnodes);
  fModel.AllNodeMap (fNodeMap);
  fModel.AllNodeIDs (nodeIDs);
  int numpoints;
  switch (selection)
    {
    case 1: // List
      {
	cout << "\n Node list defined individually\n";
	if (fWrite)
	  cout << "\n Enter the number of nodes: ";
	fIn >> numpoints;
	if (fEcho) fEchoOut << numpoints << endl;
	nodes.Dimension (numpoints);
	index.Dimension (numpoints);
	for (int n=0; n < numpoints; n++)
	  {
	    if (fWrite)
	      cout << " Enter node " << n+1 << ": ";
	    fIn >> nodes[n];
	    if (fEcho) fEchoOut << nodes[n] << endl;

	    // translate node numbers to index
	    int dex;
	    nodeIDs.HasValue (nodes[n], dex);
	    if (dex < 0 || dex >= numnodes) 
	      ExceptionT::OutOfRange ("TranslateIOManager::InitializeNodePoints","Node %i was not found", nodes[n]);
	    index [n] = dex;
	  }
	break;
      }
    case 2: // Node set
      {
	int num = fModel.NumNodeSets ();
	const ArrayT<StringT>& nodesetnames = fModel.NodeSetIDs();
	if (fWrite)
	  {
	    cout << "\n";
	    for (int h=0; h < num; h++)
	      cout << "    " << h+1 << ". " << nodesetnames[h] << "\n";
	    cout << "\n Enter the number of the node set: ";
	  }
	int ni;
	fIn >> ni;
	if (fEcho) fEchoOut << ni << endl;
	cout << "\n Node list defined by node set: " << ni << " " << nodesetnames[ni-1] << endl;
	ni--;
	numpoints = fModel.NodeSetLength (nodesetnames[ni]);
	nodes.Dimension (numpoints);
	index.Dimension (numpoints);
	index = fModel.NodeSet (nodesetnames[ni]);
	for (int n=0; n < numpoints; n++)
	  nodes[n] = nodeIDs [index[n]];
	break;
      }
    case 3: // nth node
      {
	int freq;
	if (fWrite)
	  {
	    cout << "\n Number of Nodes: " << numnodes << "\n";
	    cout << "   Enter n: ";
	  }
	fIn >> freq;
	if (fEcho) fEchoOut << freq << endl;
	cout << "\n Node list defined by every " << freq << "th node.\n";
	numpoints = numnodes/freq;
	nodes.Dimension (numpoints);
	index.Dimension (numpoints);
	for (int n=0; n < numpoints; n++)
	  {
	    nodes[n] = fNodeMap[n*freq];
	    index[n] = n*freq;
	  }
	break;
	}
	case 4: // {start, end, increment}
	{
		/* start */
		int start;
		if (fWrite) {
	    	cout << "\n Number of Nodes: " << numnodes << "\n";
	    	cout << "    Enter start node: ";
		}
		fIn >> start;
		if (fEcho) fEchoOut << start << endl;

		/* end */
		int end;
		if (fWrite) cout << "      Enter end node: ";
		fIn >> end;
		if (fEcho) fEchoOut << end << endl;

		/* increment */
		int incr;
		if (fWrite) cout << "Enter node increment: ";
		fIn >> incr;
		if (fEcho) fEchoOut << incr << endl;

		/* generate list of nodes */
		cout << "\n Node list defined by {start, end, increment} = {" 
		     << start << ",  " << end << ",  " << incr << "}\n";
		numpoints = (end - start + 1)/incr;
		nodes.Dimension (numpoints);
		index.Dimension (numpoints);
		int node = start;
		for (int n = 0; n < numpoints; n++) {
			nodes[n] = fNodeMap[node];
			index[n] = node;
			node += incr;
		}

		break;
	}
    default:
      ExceptionT::GeneralFail("TranslateIOManager::InitializeNodePoints","Invalid node list type %i", selection);
    }
}

void TranslateIOManager::SelectElements(StringT& ID, iArrayT& elements, iArrayT& index)
{
	const char caller[] = "ExtractIOManager::SelectElements";

	int selection;
	if (fWrite)
	{
		cout << "\n One file will be written per element point.\n";
		cout << "1. List of element\n";
		cout << "2. Every nth element\n";
		cout << "\n How do you want to define your list of elements: ";
	}
	fIn >> selection;
	if (fEcho) fEchoOut << selection << endl;

	int nel, nen;
	fModel.ElementGroupDimensions(ID, nel, nen);

	int num_elements;
	switch (selection)
	{
		case 1: // List
		{
			iArrayT elementIDs(nel);
			fModel.ElementIDs(ID, elementIDs);

			cout << "\n Element list defined individually\n";
			if (fWrite)
				cout << "\n Enter the number of elements: ";
			fIn >> num_elements;
			if (fEcho) fEchoOut << num_elements << endl;
			
			elements.Dimension(num_elements);
			index.Dimension(num_elements);
			elements = -1;
			for (int n=0; n < num_elements; n++)
			{
				if (fWrite) cout << " Enter element " << n+1 << ": ";
				fIn >> elements[n];
				if (fEcho) fEchoOut << elements[n] << endl;

				// translate node numbers to index
				int dex;
				elementIDs.HasValue(elements[n], dex);
				//if (dex < 0 || dex >= num_elements) 
				if (dex < 0 ) 
					ExceptionT::OutOfRange (caller, "Element %d was not found. If -1 EOF found early.", elements[n]);
				index[n] = dex;
	  		}
		break;
      	}
		case 2: // nth element
		{
			iArrayT element_map(nel);
			fModel.ElementMap(ID, element_map);

			int freq;
			if (fWrite)
			{
				cout << "\n Number of Elements: " << nel << "\n";
				cout << "   Enter n: ";
			}
			fIn >> freq;
			if (fEcho) fEchoOut << freq << endl;
			cout << "\n Element list defined by every " << freq << "th element.\n";
			int num_elem = nel/freq;
			elements.Dimension(num_elem);
			index.Dimension(num_elem);
			for (int n = 0; n < num_elem; n++)
			{
				elements[n] = element_map[n*freq];
	    		index[n] = n*freq;
			}
		break;
		}
		default:
		  {
		    cout << caller << ": Invalid selection " << selection << "\n\n";
		    ExceptionT::GeneralFail(caller);
		  }
	}
}

void TranslateIOManager::InitializeTime (void)
{
  fNumTS = fModel.NumTimeSteps ();
  fTimeSteps.Dimension (fNumTS);
  if (fNumTS > 0)
    {
      fModel.TimeSteps (fTimeSteps);

      int selection;
      if (fWrite)
	{
	  cout << "\n Number of Time Steps Available: " << fNumTS << endl;
	  if (fNumTS < 100)
	    for (int b=0; b < fNumTS; b++)
	      cout << "    " << b+1 << ". " << fTimeSteps[b] << "\n";
	  cout << "\n1. Translate All\n";
	  cout << "2. Translate Specified\n";
	  cout << "3. Translate Specified Range\n";
	  cout << "4. Translate Every nth step\n";
	  cout << "5. Translate None (just geometry)\n";
	  cout << "\n Enter Selection: ";
	}
      fIn >> selection;
      if (fEcho) fEchoOut << selection << endl;

      switch (selection)
	{
	case 1:
	  {
	    cout << "\n Translating all time steps.\n";
	    fTimeIncs.Dimension (fNumTS);
	    fTimeIncs.SetValueToPosition ();
	    break;
	  }
	case 2:
	  {
	    cout << "\n Translating list of time steps.\n";
	    if (fWrite)
	      cout << "\n Enter the number of time steps to translate: ";
	    fIn >> fNumTS;
	    if (fEcho) fEchoOut << fNumTS << endl;
	    dArrayT temp (fNumTS);
	    fTimeIncs.Dimension (fNumTS);
	    if (fWrite)
	      cout << "\n Increments are numbered consequetively from 1.\n";
	    for (int i=0; i < fNumTS; i++)
	      {
		if (fWrite)
		  cout << "    Enter time increment " << i+1 << ": ";
		fIn >> fTimeIncs[i];
		if (fEcho) fEchoOut << fTimeIncs[i] << endl;
		fTimeIncs[i]--;
		if (fTimeIncs[i] < 0 || fTimeIncs[i] >= fTimeSteps.Length())
		  throw ExceptionT::kOutOfRange;
		temp[i] = fTimeSteps[fTimeIncs[i]];
	      }
	    fTimeSteps = temp;
	    break;
	  }
	case 3:
	  {
	    int start, stop;
	    if (fWrite)
	      {
		cout << "\n Increments are numbered consequetively from 1.\n";
		cout << " Enter the starting increment: ";
	      }
	    fIn >> start;
	    if (fEcho) fEchoOut << start << endl;
	    if (fWrite)
	      cout << " Enter the end increment: ";
	    fIn >> stop;
	    if (fEcho) fEchoOut << stop << endl;
	    cout << "\n Translating time steps from " << start << " to " << stop << ".\n";
	    if (stop < start || start < 1) 
	      ExceptionT::GeneralFail ("TranslateIOManager::InitializeTime","Incorrect start/stop parameters");
	    fNumTS = stop-start+1;
	    dArrayT temp (fNumTS);
	    fTimeIncs.Dimension (fNumTS);
	    fTimeIncs.SetValueToPosition();
	    fTimeIncs += start;
	    temp.Collect (fTimeIncs, fTimeSteps);
	    fTimeSteps = temp;
	    break;
	  }
	case 4:
	  {
	    int n;
	    if (fWrite)
	      cout << "\n Enter n: ";
	    fIn >> n;
	    if (fEcho) fEchoOut << n << endl;
	    cout << "\nTranslating every " << n << "th time step.\n";
	    fNumTS = fNumTS / n;
	    fTimeIncs.Dimension (fNumTS);
	    dArrayT temp (fNumTS);
	    for (int i=0; i < fNumTS; i++)
	      {
			temp[i] = fTimeSteps [i*n];
			fTimeIncs[i] = i*n;
	      }
	    fTimeSteps = temp;
	    break;
	  }
	default:
	  {
	    fNumTS = 0;
	    fTimeSteps.Free ();
	    fTimeIncs.Free();
	  }
	}
   }
}

void TranslateIOManager::TranslateVariables(void)
{
	/* output sets */
	if (!fOutput) 
	  ExceptionT::GeneralFail ("TranslateIOManager::TranslateVariables", "Output file not initialized");

	const ArrayT<OutputSetT*>& output_sets = fOutput->ElementSets();
	if (output_sets.Length() != fOutputID.Length()) throw ExceptionT::kSizeMismatch;

	const ArrayT<StringT>& names = fModel.ElementGroupIDs ();	

	/* loop over output groups */
	for (int g = 0; g < output_sets.Length(); g++)
	{
		OutputSetT& set = *(output_sets[g]);
		int num_nodes = set.NumNodes();
		int num_node_values = set.NumNodeValues();
		int num_elems = set.NumElements();
		int num_elem_values = set.NumElementValues();

		/* work space */
		dArray2DT n_values(num_nodes, fNumNV);
		dArray2DT e_values(num_elems, fNumEV);

		/* filtering nodal values */
		dArray2DT n_values_all;
		if (fModel.NumNodeVariables() != fNumNV)
			n_values_all.Dimension(num_nodes, fModel.NumNodeVariables());
		else
			n_values_all.Alias(n_values);

		/* filtering element values */
		dArray2DT e_values_all;
		if (fModel.NumElementVariables() != fNumEV)
			e_values_all.Dimension(num_elems, fModel.NumElementVariables());
		else
			e_values_all.Alias(e_values);
	
		/* loop over time steps */
		for (int t = 0; t < fTimeIncs.Length(); t++)
		{
			/* report */
			if (g == 0) cout << "Time Step " << fTimeIncs[t]+1 << ": " << fTimeSteps[t] << endl;

			if (fOneOutputSet)
			  {
			    /* read node values */
			    fModel.AllNodeVariables(fTimeIncs[t], n_values_all);
			    
			    /* read values for all blocks - assumes block values assembled
			     * one after the next */
			    fModel.AllElementVariables(fTimeIncs[t], e_values_all);
			  }
			else
			  {
			    /* read node values */
			    fModel.NodeVariables (fTimeIncs[t], names[g], n_values_all);

			    /* read element values */
			    fModel.ElementVariables (fTimeIncs[t], names[g], e_values_all);
			  }
			  
			/* filter nodal values */
			if (n_values_all.MinorDim() != n_values.MinorDim()) {
				int index = 0;
				for (int i = 0; i < n_values_all.MinorDim(); i++)
					if (fKeepNodeLabels[i]) 
						n_values.ColumnCopy(index++, n_values_all, i);	
			}

			/* filter nodal values */
			if (e_values_all.MinorDim() != e_values.MinorDim()) {
				int index = 0;
				for (int i = 0; i < e_values_all.MinorDim(); i++)
					if (fKeepElementLabels[i]) 
						e_values.ColumnCopy(index++, e_values_all, i);	
			}

			/* write it */
			fOutput->WriteOutput(fTimeSteps[t], fOutputID[g], n_values, e_values);
		}
	}
}

void TranslateIOManager::WriteGeometry (void)
{
  WriteNodes ();
  WriteNodeSets ();
  WriteElements ();
  WriteSideSets ();
}

void TranslateIOManager::WriteNodes (void)
{
	int nnd = fModel.NumNodes();
	int nsd = fModel.NumDimensions();
	fNodeMap.Dimension(nnd);
	fModel.AllNodeMap(fNodeMap);
	fModel.AllNodeIDs(fNodeID);
	/*fNodeMap--; offset done by model manager */
	/* Difference between "node map" and "nodes used". 
	 * 
	 * nodes used refers to nodes used by an element group
	 * node map is just an global index of nodes in fCoordinates
	 *
	 */
  
	/* node map should not be empty */
	if (fNodeMap.Length() == 0)
	  ExceptionT::GeneralFail ("TranslateIOManager::WriteNodes","Node number map is empty.");
  
	/* do not need to do any mapping, fNodeMap is global and offset 
	   when it comes from modelmanager */
  
	fOutput->SetCoordinates (fModel.Coordinates(), &fNodeID);
	cout << "\n Number of Nodes: " << nnd << " dim: " << nsd << endl;
}

void TranslateIOManager::WriteNodeSets (void)
{
  int num = fModel.NumNodeSets ();
  if (num <= 0) return;
  const ArrayT<StringT>& names = fModel.NodeSetIDs();

  int selection;
  cout << "\n Number of Node Sets: " << num << endl;
  if (fWrite)
    {
      cout << "\n1. Translate All\n";
      cout << "2. Translate Some\n";
      cout << "3. Translate None\n";
      cout << "\n selection: ";
    }
  fIn >> selection;
  if (fEcho) fEchoOut << selection << endl;

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer("yes");
      if (selection == 2)
	{
	  if (fWrite)
	    cout << "    Translate Node Set " << names[i] << " (y/n) ? ";
	  fIn >> answer;
	  if (fEcho) fEchoOut << answer[0] << endl;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y') 
	fOutput->AddNodeSet (fModel.NodeSet(names[i]), names[i]);
    }
}
 
void TranslateIOManager::WriteElements(void)
{
  // number of sets
	int num = fModel.NumElementGroups ();
	if (num == 0) 
	  ExceptionT::GeneralFail ("TranslateIOManager::WriteElements","No element sets");

	// which to output
	cout << "\n Number of Element Groups: " << num << endl;
	if (fWrite)
	  {
	    cout << "\n1. Translate All\n";
	    cout << "2. Translate Some\n";
	    cout << "\n selection: ";
	  }
	int selection;
	fIn >> selection;
	if (fEcho) fEchoOut << selection << endl;

	/* work space */
	AutoArrayT<const iArray2DT*> blocks;
	AutoArrayT<StringT> block_IDs;
	bool changing = false;
	StringT answer;
	const ArrayT<StringT>& names = fModel.ElementGroupIDs();
	GeometryT::CodeT geometry = GeometryT::kNone;

	// see if all element groups have the same geometry code
	// some databases have element sets of differing geometry codes
	fOneOutputSet = true;
	geometry = fModel.ElementGroupGeometry(names[0]);
	for (int h = 1; h < num && fOneOutputSet; h++)
	  if (geometry != fModel.ElementGroupGeometry (names[h]))
	    fOneOutputSet = false;
	if (fOneOutputSet)
	  fOutputID.Dimension (1);
	else
	  fOutputID.Dimension (num);

	for (int e = 0; e < num; e++)
	  {
		// allow user to select element groups
		answer = "yes";
		if (selection == 2) 
		  {
		    if (fWrite) cout << "    Translate Element Group " << names[e] << " (y/n) ? ";
		    fIn >> answer;
		    if (fEcho) fEchoOut << answer[0] << endl;
		  }
      
		if (answer [0] == 'y' || answer[0] == 'Y')
		  {
			/* read connectivities */
			fModel.ReadConnectivity(names[e]);

			if (fOneOutputSet)
			  {
			    // if same geocode, save data for sending to outputset later
			    blocks.Append(fModel.ElementGroupPointer(names[e]));
			    block_IDs.Append(names[e]);
			  }
			else if (fModel.ElementGroup(names[e]).MajorDim() > 0) /* skip empty sets */
			  {
			    // if differing codes, send output now
			    ArrayT<StringT> block_ID(1);
			    ArrayT<const iArray2DT*> conn (1);
			    block_ID[0] = names[e];
			    StringT ID;
			    ID.Append(e+1);
			    conn[0] = fModel.ElementGroupPointer (names[e]);
			    OutputSetT set(fModel.ElementGroupGeometry (names[e]), block_ID, 
					   conn, fNodeLabels, fElementLabels, changing);

			    fOutputID[e] = fOutput->AddElementSet (set);
			  }
			else
			  cout << " Element Group " << names[e] << " has no elements.\n";
		  }
	  }

	if (fOneOutputSet)
	  {
	    /* set ID */
	    StringT ID;
	    if (block_IDs.Length() == 1)
	      ID = block_IDs[0];
	    else
	      ID.Append(0);
	    
	    /* create output set */
	    OutputSetT set(geometry, block_IDs, blocks, fNodeLabels, fElementLabels, changing);
	    fOutputID[0] = fOutput->AddElementSet(set);
	  }
}

void TranslateIOManager::WriteSideSets (void)
{
  int num = fModel.NumSideSets ();
  if (num <= 0) return;
  fGlobalSideSets.Dimension (num);
  const ArrayT<StringT>& names = fModel.SideSetIDs();

  int selection;
  cout << "\n Number of Side Sets: " << num << endl;
  if (fWrite)
    {
      cout << "\n1. Translate All\n";
      cout << "2. Translate Some\n";
      cout << "3. Translate None\n";
      cout << "\n selection: ";
    }
  fIn >> selection;
  if (fEcho) fEchoOut << selection << endl;

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer ("yes");
      if (selection == 2)
	{
	  if (fWrite)
	    cout << "    Translate Side Set " << names[i] << " (y/n) ? ";
	  fIn >> answer;
	  if (fEcho) fEchoOut << answer[0] << endl;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y')
	{
	  fGlobalSideSets[i] = fModel.SideSet (names[i]);
	  const StringT& g = fModel.SideSetGroupID(names[i]);
	    
	  fOutput->AddSideSet (fGlobalSideSets [i], names[i], g);
	}
    }
}

void TranslateIOManager::VariableQuery (const ArrayT<StringT>& names, iArrayT& list)
{
  iAutoArrayT temp;
  for (int i=0; i < names.Length(); i++)
    {
      StringT answer;
      if (fWrite)
	cout << " Extract variable " << names[i] << " (y/n) ? ";
      fIn >> answer;
      if (fEcho) fEchoOut << answer[0] << endl;

      if (answer[0] == 'y' || answer[0] == 'Y')
	temp.Append (i);
    }

  list.Dimension (temp.Length());
  list.CopyPart (0, temp, 0, temp.Length());  
}

void TranslateIOManager::ReNameLabels(const StringT& data_type, ArrayT<StringT>& labels)
{
	if (labels.Length() == 0) return;

	if (fWrite)
	  cout << "\n Rename " << data_type << " labels (y/n) ? ";
	StringT reply;
	fIn >> reply;
	if (fEcho) fEchoOut << reply[0] << endl;

	/* clear newline */
	fstreamT::ClearLine(fIn);

	if (reply[0] == 'y' || reply[0] == 'Y')
		for (int i = 0; i < labels.Length(); i++)
		{
		  if (fWrite)
			cout << " Rename " << labels[i] << " <" << labels[i] << "> : ";

			/* peek at reply */
			char test = fIn.peek();

			/* not empty */
			if (test != '\n') {
				fIn >> reply;
				labels[i] = reply;
				if (fEcho) fEchoOut << reply << endl;
			} 	
			
			/* clear line */
			fstreamT::ClearLine(fIn);
		}
}

void TranslateIOManager::FilterLabels(const StringT& data_type, iArrayT& filter,
	ArrayT<StringT>& labels)
{
	if (labels.Length() == 0) return;

	if (fWrite) cout << "\n Filter " << data_type << " labels (y/n) ? ";
	StringT reply;
	fIn >> reply;
	if (fEcho) fEchoOut << reply[0] << endl;

	/* clear newline */
	fstreamT::ClearLine(fIn);

	if (reply[0] == 'y' || reply[0] == 'Y')
	{
		filter.Dimension(labels.Length());
		filter = 1;
		for (int i = 0; i < labels.Length(); i++)
		{
		  if (fWrite)
			cout << " Keep " << labels[i] << " <y> ? ";

			/* peek at reply */
			char test = fIn.peek();

			/* not empty */
			if (test != '\n') {
				fIn >> reply;
				if (reply[0] == 'n' || reply[0] == 'N')
					filter[i] = 0;
				if (fEcho) fEchoOut << reply << endl;
			}
			
			/* clear line */
			fstreamT::ClearLine(fIn);
		}
		
		/*  filter labels */
		int keep_count = filter.Count(1);
		ArrayT<StringT> keep_labels(keep_count);
		keep_count = 0;
		for (int i = 0; i < filter.Length(); i++)
			if (filter[i])
				keep_labels[keep_count++] = labels[i];
				
		/* return filtered list */
		labels.Swap(keep_labels);
	}
}
