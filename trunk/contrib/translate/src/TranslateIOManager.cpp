/* $Id: TranslateIOManager.cpp,v 1.15 2002-02-18 09:44:06 paklein Exp $  */

#include "TranslateIOManager.h"
#include "IOBaseT.h"
#include "OutputSetT.h"

#include "AbaqusOutputT.h"
#include "AVSOutputT.h"
#include "EnSightOutputT.h"
#include "ExodusOutputT.h"
#include "TecPlotOutputT.h"
#include "FE_ASCIIT.h"

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
  fCoords (0)
{
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
    {
      cout << "\n Unable to initialize model file\n";
      throw eGeneralFail;
    }
}

void TranslateIOManager::SetOutput (const StringT& program_name, const StringT& version, const StringT& title)
{
  IOBaseT temp (cout);
  int outputformat = -1;
  if (fWrite)
    {
      cout << "\n\n";
      temp.OutputFormats (cout);
      cout << "\n Enter the Output Format: ";
    }
  fIn >> outputformat;

  if (fWrite)
    cout << "\n Enter the root of the output files: ";
  fIn >> fOutputName;
  cout << "\n Output format: " << outputformat << " File: " << fOutputName << endl;

  fOutputName.ToNativePathName();
  fOutputName.Append(".ext"); //trimmed off by fOutput

  ArrayT<StringT> outstrings (4);
  outstrings[0] = fOutputName;
  outstrings[1] = title;
  outstrings[2] = program_name;
  outstrings[3] = version;

  switch (outputformat)
    {
    case IOBaseT::kExodusII:
      fOutput = new ExodusOutputT (fMessage, outstrings);
      break;
    case IOBaseT::kTahoeII:
      fOutput = new FE_ASCIIT (fMessage, true, outstrings);
      break;
    case IOBaseT::kEnSight:
      fOutput = new EnSightOutputT (fMessage, outstrings, 4, false);
      break;
    case IOBaseT::kEnSightBinary:
      fOutput = new EnSightOutputT (fMessage, outstrings, 4, true);
      break;
    case IOBaseT::kAbaqus:
      fOutput = new AbaqusOutputT (fMessage, outstrings, false);
      break;
    case IOBaseT::kAbaqusBinary:
      fOutput = new AbaqusOutputT (fMessage, outstrings, true);
      break;
    case IOBaseT::kTecPlot:
      fOutput = new TecPlotOutputT (fMessage, outstrings, 4);
      break;
    case IOBaseT::kAVS:
    case IOBaseT::kAVSBinary:
      fOutput = new AVSOutputT (fMessage, outstrings, false);
      break;
    default:
      {
	fMessage << "\n Unknown output format: " << outputformat << "\n";
	throw eDatabaseFail;
      }
    }
  if (!fOutput) throw eOutOfMemory;
}

void TranslateIOManager::InitializeVariables (void)
{
  fNumNV = fModel.NumNodeVariables ();
  fNumEV = fModel.NumElementVariables ();
  fNumQV = fModel.NumQuadratureVariables ();

  cout << "\n" << setw (10) << fNumNV << " Node Variables\n";
  cout << setw (10) << fNumEV << " Element Variables\n";
  cout << setw (10) << fNumQV << " Quadrature Varaibles\n";

  fNodeLabels.Allocate (fNumNV);
  fElementLabels.Allocate (fNumEV);

  if (fNumNV > 0) fModel.NodeLabels (fNodeLabels);
  if (fNumEV > 0) fModel.ElementLabels (fElementLabels);

  // future: query user as to which variables to translate
}

void TranslateIOManager::InitializeNodeVariables (void)
{
  fNumNV = fModel.NumNodeVariables ();
  cout << "\n" << setw (10) << fNumNV << " Node Variables\n\n";
  fNodeLabels.Allocate (fNumNV);
  if (fNumNV > 0) fModel.NodeLabels (fNodeLabels);

  // query user as to which variables to translate
  VariableQuery (fNodeLabels, fNVUsed);

  StringT answer;
  if (fWrite)
    cout << "\n Do you wish to translate coordinate values (y/n) ? ";
  fIn >> answer;
  
  if (answer[0] == 'y' || answer[0] == 'Y')
    {
      fCoords = true;
      int numnodes;
      fModel.CoordinateDimensions (numnodes, fCoords);
      cout << "\n Adding coordinates to variable values.\n";
    }
  else
    fCoords = 0;
}

void TranslateIOManager::InitializeQuadVariables (void)
{
  fNumQV = fModel.NumQuadratureVariables ();
  cout << "\n" << setw (10) << fNumQV << " Quadrature Variables\n\n";
  fQuadratureLabels.Allocate (fNumQV);
  if (fNumQV > 0) fModel.QuadratureLabels (fQuadratureLabels);

  // query user as to which variables to translate
  if (fNumQV < 1)
    {
      fMessage << "\n No quadrature variables found.";
      throw eGeneralFail;
    }
  //cout << fNumQV << " " << fQuadratureLabels[0] <<endl;
  VariableQuery (fQuadratureLabels, fQVUsed);
}

void TranslateIOManager::InitializeElements (int& group, StringT& groupname) const
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
  if (group < 1 || group > elemsetnames.Length()) 
    {
      cout << "\n The number entered for an element group is invalid: "
	   << group << "\n";
      cout << "Minimum limit is 1 and maximum is " << elemsetnames.Length() << endl;
      throw eOutOfRange;
    }
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
      cout << "\n How do you want to define your list of nodes: ";
    }
  fIn >> selection;

  int numnodes, numdims;
  fModel.CoordinateDimensions (numnodes, numdims);
  fNodeMap.Allocate (numnodes);
  fModel.AllNodeMap (fNodeMap);
  int numpoints;
  switch (selection)
    {
    case 1:
      {
	cout << "\n Node list defined individually\n";
	if (fWrite)
	  cout << "\n Enter the number of nodes: ";
	fIn >> numpoints;
	nodes.Allocate (numpoints);
	index.Allocate (numpoints);
	for (int n=0; n < numpoints; n++)
	  {
	    if (fWrite)
	      cout << " Enter node " << n+1 << ": ";
	    fIn >> nodes[n];

	    // translate node numbers to index
	    int dex;
	    fNodeMap.HasValue (nodes[n], dex);
	    if (dex < 0 || dex >= numnodes) 
	      {
		cout << " ExtractIOManager::InitializeNodePoints\n";
		cout << " Node " << nodes[n] << " was not found.\n";
		throw eOutOfRange;
	      }
	    index [n] = dex;
	  }
	break;
      }
    case 2:
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
	cout << "\n Node list defined by node set: " << ni << " " << nodesetnames[ni-1] << endl;
	ni--;
	numpoints = fModel.NodeSetLength (nodesetnames[ni]);
	nodes.Allocate (numpoints);
	index.Allocate (numpoints);
	index = fModel.NodeSet (nodesetnames[ni]);
	for (int n=0; n < numpoints; n++)
	  nodes[n] = fNodeMap [index[n]];
	break;
      }
    case 3:
      {
	int freq;
	if (fWrite)
	  {
	    cout << "\n Number of Nodes: " << numnodes << "\n";
	    cout << "   Enter n: ";
	  }
	fIn >> freq;
	cout << "\n Node list defined by every " << freq << "th node.\n";
	numpoints = numnodes/freq;
	nodes.Allocate (numpoints);
	index.Allocate (numpoints);
	for (int n=0; n < numpoints; n++)
	  {
	    nodes[n] = fNodeMap[n*freq];
	    index[n] = n*freq;
	  }
	break;
      }
    default:
      throw eGeneralFail;
    }
}

void TranslateIOManager::InitializeTime (void)
{
  fNumTS = fModel.NumTimeSteps ();
  fTimeSteps.Allocate (fNumTS);
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

      switch (selection)
	{
	case 1:
	  {
	    cout << "\n Translating all time steps.\n";
	    fTimeIncs.Allocate (fNumTS);
	    fTimeIncs.SetValueToPosition ();
	    break;
	  }
	case 2:
	  {
	    cout << "\n Translating list of time steps.\n";
	    if (fWrite)
	      cout << "\n Enter the number of time steps to translate: ";
	    fIn >> fNumTS;
	    dArrayT temp (fNumTS);
	    fTimeIncs.Allocate (fNumTS);
	    if (fWrite)
	      cout << "\n Increments are numbered consequetively from 1.\n";
	    for (int i=0; i < fNumTS; i++)
	      {
		if (fWrite)
		  cout << "    Enter time increment " << i+1 << ": ";
		fIn >> fTimeIncs[i];
		fTimeIncs[i]--;
		if (fTimeIncs[i] < 0 || fTimeIncs[i] >= fTimeSteps.Length())
		  throw eOutOfRange;
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
	    if (fWrite)
	      cout << " Enter the end increment: ";
	    fIn >> stop;
	    cout << "\n Translating time steps from " << start << " to " << stop << ".\n";
	    if (stop < start) throw eGeneralFail;
	    if (start < 1) throw eGeneralFail;
	    fNumTS = stop-start+1;
	    dArrayT temp (fNumTS);
	    fTimeIncs.Allocate (fNumTS);
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
	    cout << "\nTranslating every " << n << "th time step.\n";
	    fNumTS = fNumTS / n;
	    fTimeIncs.Allocate (fNumTS);
	    dArrayT temp (fNumTS);
	    for (int i=0; i < fNumTS; i++)
	      {
		temp[i] = fTimeSteps [i*n];
		fTimeIncs = i*n;
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
	if (!fOutput) {
		cout << "\n TranslateIOManager::TranslateVariables: output not initialized" << endl;
		throw eGeneralFail;
	}
	const ArrayT<OutputSetT*>& output_sets = fOutput->ElementSets();
	if (output_sets.Length() != fOutputID.Length()) throw eSizeMismatch;

	/* loop over element groups */
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
	
		/* loop over time steps */
		for (int t = 0; t < fTimeIncs.Length(); t++)
		{
			/* report */
			if (g == 0) cout << "Time Step " << fTimeIncs[t]+1 << ": " << fTimeSteps[t] << endl;

			/* read node values */
			fModel.AllNodeVariables(t, n_values);
		
			/* read values for all blocks - assumes block values assembled
			 * one after the next */
			fModel.AllElementVariables(t, e_values);

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
	int nnd, nsd;
	fModel.CoordinateDimensions(nnd, nsd);
	fNodeMap.Dimension(nnd);
	fModel.AllNodeMap(fNodeMap);
	fNodeMap--;
  
	/* node map should not be empty */
	if (fNodeMap.Length() == 0) {
		cout << "\n TranslateIOManager::WriteNodes: node number map is empty" << endl;
		throw eGeneralFail;
		}
  
	/* Difference between "node map" and "nodes used". Coordinates list is "global" if:
	 *
     *   fCoordinates[fNodeMap[i]] == coordinates of node fNodeMap[i] 
     *
     */
	int min, max;
	fNodeMap.MinMax(min, max);
	if (min != 0 || max >= nnd) /* this is not a global node set */
	{
		fCoordinates.Dimension(max + 1, nsd);
		fCoordinates.Assemble(fNodeMap, fModel.Coordinates());
	}
	else fCoordinates.Alias(fModel.Coordinates());
  
	fOutput->SetCoordinates (fCoordinates, &fNodeMap);
	cout << "\n Number of Nodes: " << nnd << " dim: " << nsd << endl;
}

void TranslateIOManager::WriteNodeSets (void)
{
  int num = fModel.NumNodeSets ();
  if (num <= 0) return;
  const ArrayT<StringT>& names = fModel.NodeSetIDs();

  int selection;
  if (fWrite)
    {
      cout << "\n Number of Node Sets: " << num << endl;
      cout << "\n1. Translate All\n";
      cout << "2. Translate Some\n";
      cout << "3. Translate None\n";
      cout << "\n selection: ";
    }
  fIn >> selection;

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer("yes");
      if (selection == 2)
	{
	  if (fWrite)
	    cout << "    Translate Node Set " << names[i] << " (y/n) ? ";
	  fIn >> answer;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y') 
	{
	  int ID = atoi(names[i]);
	  fOutput->AddNodeSet (fModel.NodeSet(names[i]), ID);
	}
    }
}
 
void TranslateIOManager::WriteElements(void)
{
	int num = fModel.NumElementGroups ();
	cout << "\n Number of Element Groups: " << num << endl;
	cout << "\n1. Translate All\n";
	cout << "2. Translate Some\n";
	cout << "\n selection: ";
	int selection;
	cin >> selection;

	/* element sets */
	AutoArrayT<const iArray2DT*> blocks;
	AutoArrayT<StringT> block_IDs;

	bool changing = false;
	StringT answer;
	const ArrayT<StringT>& names = fModel.ElementGroupIDs();
	if (num == 0) {
		cout << "\n TranslateIOManager::WriteElements: no element sets" << endl;
		throw eGeneralFail;
	}
	GeometryT::CodeT geometry = GeometryT::kNone;
	for (int e = 0; e < num; e++)
	{
		// allow user to select element groups
		answer = "yes";
		if (selection == 2) {
			cout << "    Translate Element Group " << names[e] << " (y/n) ? ";
			cin >> answer;
		}
      
		if (answer [0] == 'y' || answer[0] == 'Y')
		{
			/* read connectivities */
			fModel.ReadConnectivity(names[e]);
			blocks.Append(fModel.ElementGroupPointer(names[e]));
			block_IDs.Append(names[e]);

			/* must all be the same geometry */
			if (fModel.ElementGroup(names[e]).MajorDim() > 0) /* skip empty sets */
			{
				if (geometry == GeometryT::kNone)
					geometry = fModel.ElementGroupGeometry(names[e]);
				else if (geometry != fModel.ElementGroupGeometry(names[e])) {
					cout << "\n TranslateIOManager::WriteElements: mismatched geometry codes\n" 
					     <<   "     " << geometry <<  " and " << fModel.ElementGroupGeometry(names[e]) 
					     << endl;
					throw eDatabaseFail;			
				}
			}

#if 0
			if (true || conn[0]->Length() > 0)
			{
				ArrayT<StringT> block_ID(1);
				block_ID[0] = names[e];
				StringT ID;
				ID.Append(e+1);
//				OutputSetT set(ID, fModel.ElementGroupGeometry (names[e]), block_ID, 
//					conn, fNodeLabels, fElementLabels, changing);
//				fOutputID[e] = fOutput->AddElementSet (set);
			}
			else
				cout << " Element Group " << names[e] << " has no elements.\n";
#endif
		}
	}
	
	/* set ID */
	StringT ID;
	if (block_IDs.Length() == 1)
		ID = block_IDs[0];
	else
		ID.Append(0);

	/* create output set */
	fOutputID.Allocate(1);
	OutputSetT set(ID, geometry, block_IDs, blocks, fNodeLabels, fElementLabels, changing);
	fOutputID[0] = fOutput->AddElementSet(set);
}

void TranslateIOManager::WriteSideSets (void)
{
  int num = fModel.NumSideSets ();
  if (num <= 0) return;
  fGlobalSideSets.Allocate (num);
  const ArrayT<StringT>& names = fModel.SideSetIDs();

  int selection;
  if (fWrite)
    {
      cout << "\n Number of Side Sets: " << num << endl;
      cout << "\n1. Translate All\n";
      cout << "2. Translate Some\n";
      cout << "3. Translate None\n";
      cout << "\n selection: ";
    }
  fIn >> selection;

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer ("yes");
      if (selection == 2)
	{
	  if (fWrite)
	    cout << "    Translate Node Set " << names[i] << " (y/n) ? ";
	  fIn >> answer;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y')
	{
	  fGlobalSideSets[i] = fModel.SideSet (names[i]);
	  const StringT& g = fModel.SideSetGroupID(names[i]);
	    
	  int g_int = atoi(g);
	  int ID = atoi(names[i]);
	  fOutput->AddSideSet (fGlobalSideSets [i], ID, g_int);
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

      if (answer[0] == 'y' || answer[0] == 'Y')
	temp.Append (i);
    }

  list.Allocate (temp.Length());
  list.CopyPart (0, temp, 0, temp.Length());  
}
