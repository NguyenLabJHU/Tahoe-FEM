/* $Id: TranslateIOManager.cpp,v 1.31 2002-10-28 14:19:02 sawimme Exp $  */
#include "TranslateIOManager.h"

#include "ExceptionT.h"
#include "IOBaseT.h"
#include "OutputSetT.h"
#include "AbaqusOutputT.h"
#include "AVSOutputT.h"
#include "EnSightOutputT.h"
#include "ExodusOutputT.h"
#include "TecPlotOutputT.h"
#include "FE_ASCIIT.h"
#include "PatranOutputT.h"

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
      throw ExceptionT::kGeneralFail;
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
  cout << "\n Output format: " << outputformat << "\n File: " << fOutputName << endl;

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
      case IOBaseT::kPatranNeutral:
	fOutput = new PatranOutputT (fMessage, outstrings, false);
	break;
    case IOBaseT::kAVS:
    case IOBaseT::kAVSBinary:
      fOutput = new AVSOutputT (fMessage, outstrings, false);
      break;
    default:
      {
	fMessage << "\n Unknown output format: " << outputformat << "\n";
	throw ExceptionT::kDatabaseFail;
      }
    }
  if (!fOutput) throw ExceptionT::kOutOfMemory;
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
  fElementLabels.Dimension (fNumEV);

  if (fNumNV > 0) {
  	fModel.NodeLabels (fNodeLabels);
  	ReNameLabels("node", fNodeLabels);
  }
  
  if (fNumEV > 0) {
  	fModel.ElementLabels (fElementLabels);
  	ReNameLabels("element", fElementLabels);
	}

  // future: query user as to which variables to translate
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
  
  if (answer[0] == 'y' || answer[0] == 'Y')
    {
      fCoords = fModel.NumDimensions();
      int numnodes = fModel.NumNodes();
      cout << "\n Adding coordinates to variable values.\n";
    }
  else
    fCoords = 0;
}

void TranslateIOManager::InitializeQuadVariables (void)
{
  fNumQV = fModel.NumQuadratureVariables ();
  cout << "\n" << setw (10) << fNumQV << " Quadrature Variables\n\n";
  fQuadratureLabels.Dimension (fNumQV);
  if (fNumQV > 0) fModel.QuadratureLabels (fQuadratureLabels);

  // query user as to which variables to translate
  if (fNumQV < 1)
    {
      fMessage << "\n No quadrature variables found.";
      throw ExceptionT::kGeneralFail;
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
      throw ExceptionT::kOutOfRange;
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
	nodes.Dimension (numpoints);
	index.Dimension (numpoints);
	for (int n=0; n < numpoints; n++)
	  {
	    if (fWrite)
	      cout << " Enter node " << n+1 << ": ";
	    fIn >> nodes[n];

	    // translate node numbers to index
	    int dex;
	    nodeIDs.HasValue (nodes[n], dex);
	    if (dex < 0 || dex >= numnodes) 
	      {
		cout << " ExtractIOManager::InitializeNodePoints\n";
		cout << " Node " << nodes[n] << " was not found.\n";
		throw ExceptionT::kOutOfRange;
	      }
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
    default:
      throw ExceptionT::kGeneralFail;
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
	    dArrayT temp (fNumTS);
	    fTimeIncs.Dimension (fNumTS);
	    if (fWrite)
	      cout << "\n Increments are numbered consequetively from 1.\n";
	    for (int i=0; i < fNumTS; i++)
	      {
		if (fWrite)
		  cout << "    Enter time increment " << i+1 << ": ";
		fIn >> fTimeIncs[i];
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
	    if (fWrite)
	      cout << " Enter the end increment: ";
	    fIn >> stop;
	    cout << "\n Translating time steps from " << start << " to " << stop << ".\n";
	    if (stop < start) throw ExceptionT::kGeneralFail;
	    if (start < 1) throw ExceptionT::kGeneralFail;
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
	    cout << "\nTranslating every " << n << "th time step.\n";
	    fNumTS = fNumTS / n;
	    fTimeIncs.Dimension (fNumTS);
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
	if (!fOutput) 
	  {
	    cout << "\n TranslateIOManager::TranslateVariables: output not initialized" << endl;
	    throw ExceptionT::kGeneralFail;
	  }

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
	
		/* loop over time steps */
		for (int t = 0; t < fTimeIncs.Length(); t++)
		{
			/* report */
			if (g == 0) cout << "Time Step " << fTimeIncs[t]+1 << ": " << fTimeSteps[t] << endl;

			if (fOneOutputSet)
			  {
			    /* read node values */
			    fModel.AllNodeVariables(fTimeIncs[t], n_values);
			    
			    /* read values for all blocks - assumes block values assembled
			     * one after the next */
			    fModel.AllElementVariables(fTimeIncs[t], e_values);
			  }
			else
			  {
			    /* read node values */
			    fModel.NodeVariables (fTimeIncs[t], names[g], n_values);

			    /* read element values */
			    fModel.ElementVariables (fTimeIncs[t], names[g], e_values);
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
	if (fNodeMap.Length() == 0) {
		cout << "\n TranslateIOManager::WriteNodes: node number map is empty" << endl;
		throw ExceptionT::kGeneralFail;
		}
  
	/* do not need to do any mapping, fNodeMap is global and offset 
	   when it comes from modelmanager */
	/*int min, max;
	  fNodeMap.MinMax(min, max);
	  if (min != 0 || max >= nnd)
	  {
	  fCoordinates.Dimension(max + 1, nsd);
	  fCoordinates.Assemble(fNodeMap, fModel.Coordinates());
	  fCoordinates.WriteNumbered (cout);
	  throw ExceptionT::kGeneralFail;
	  }
	  else fCoordinates.Alias(fModel.Coordinates()); */
  
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
	fOutput->AddNodeSet (fModel.NodeSet(names[i]), names[i]);
    }
}
 
void TranslateIOManager::WriteElements(void)
{
  // number of sets
	int num = fModel.NumElementGroups ();
	if (num == 0) 
	  {
	    cout << "\n TranslateIOManager::WriteElements: no element sets" << endl;
	    throw ExceptionT::kGeneralFail;
	  }

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

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer ("yes");
      if (selection == 2)
	{
	  if (fWrite)
	    cout << "    Translate Side Set " << names[i] << " (y/n) ? ";
	  fIn >> answer;
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

	/* clear newline */
	char line[255];
	fIn.getline(line, 254);

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
			} 	
			
			/* clear line */
			fIn.getline(line, 254);
		}
}
