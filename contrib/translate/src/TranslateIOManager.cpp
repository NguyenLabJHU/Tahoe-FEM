/* $Id: TranslateIOManager.cpp,v 1.5 2001-09-21 15:49:25 sawimme Exp $  */

#include "TranslateIOManager.h"
#include "IOBaseT.h"
#include "OutputSetT.h"

#include "AbaqusOutputT.h"
#include "AVSOutputT.h"
#include "EnSightOutputT.h"
#include "ExodusOutputT.h"
#include "TecPlotOutputT.h"
#include "FE_ASCIIT.h"

TranslateIOManager::TranslateIOManager (ostream& out) :
  fMessage (out),
  fModel (out),
  fOutput (NULL),
  fNumNV (0),
  fNumEV (0),
  fNumQV (0),
  fNumTS (0)
{
}

void TranslateIOManager::Translate (const StringT& program, const StringT& version, const StringT& title)
{
  fModel.Initialize ();
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
      fMessage << "\n No time steps found, writing geometry only.";
      WriteGeometry ();
      fOutput->WriteGeometry ();
      return;      
    }
  
  TranslateVariables ();
}


/**************** PRIVATE **********************/

void TranslateIOManager::SetOutput (const StringT& program_name, const StringT& version, const StringT& title)
{
  IOBaseT temp (cout);
  cout << "\n\n";
  temp.OutputFormats (cout);
  int outputformat;
  cout << "\n Enter the Output Format: ";
  cin >> outputformat;
  cout << "\n Enter the root of the output files: ";
  cin >> fOutputName;

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

void TranslateIOManager::InitializeTime (void)
{
  fNumTS = fModel.NumTimeSteps ();
  fTimeSteps.Allocate (fNumTS);
  if (fNumTS > 0)
    {
      fModel.TimeSteps (fTimeSteps);

      int selection;
      cout << "\n Number of Time Steps Available: " << fNumTS << endl;
      if (fNumTS < 100)
	for (int b=0; b < fNumTS; b++)
	  cout << "    " << b+1 << ". " << fTimeSteps[b] << "\n";
      cout << "\n1. Translate All\n";
      cout << "2. Translate Specified\n";
      cout << "3. Translate Every nth step\n";
      cout << "4. Translate None (just geometry)\n";
      cout << "\n Enter Selection: ";
      cin >> selection;

      switch (selection)
	{
	case 1:
	  {
	    fTimeIncs.Allocate (fNumTS);
	    fTimeIncs.SetValueToPosition ();
	    break;
	  }
	case 2:
	  {
	    cout << "\n Enter the number of time steps to translate: ";
	    cin >> fNumTS;
	    dArrayT temp (fNumTS);
	    fTimeIncs.Allocate (fNumTS);
	    cout << "\n Increments are numbered consequetively from 1.\n";
	    for (int i=0; i < fNumTS; i++)
	      {
		cout << "    Enter time increment " << i+1 << ": ";
		cin >> fTimeIncs[i];
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
	    int n;
	    cout << "\n Enter n: ";
	    cin >> n;
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

void TranslateIOManager::TranslateVariables (void)
{
  int ng = fModel.NumElementGroups ();
  if (ng < 1)
    {
      fMessage << "\n No element groups found.\n";
      return;
    }
  else if (ng != fOutputID.Length()) 
    throw eSizeMismatch;

  ArrayT<StringT> groupnames (ng);
  fModel.ElementGroupNames (groupnames);

  int numnodes, numdims;
  fModel.CoordinateDimensions (numnodes, numdims);
  for (int t=0; t < fNumTS; t++)
    {
      cout << "Time Step " << fTimeIncs[t]+1 << ": " << fTimeSteps[t] << endl;
      for (int g=0; g < fOutputID.Length(); g++)
	{
	  // account for user specified element groups not to translate
	  if (fOutputID [g] > -1) 
	    {
	      int numelems, numelemnodes;
	      fModel.ElementGroupDimensions (g, numelems, numelemnodes);
	      
	      // read all variables for this step
	      dArray2DT nvalues (numnodes, fNumNV); // numnodes is larger than needed
	      dArray2DT evalues (numelems, fNumEV);
	      if (fNumNV > 0) 
		fModel.NodeVariables (fTimeIncs[t], groupnames[g], nvalues);
	      if (fNumEV > 0) 
		fModel.ElementVariables (fTimeIncs[t], groupnames[g], evalues);
	      
	      // future: accout for which variables user wanted translated
	      
	      fOutput->WriteOutput (fTimeSteps[t], fOutputID[g], nvalues, evalues);
	    }
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
  int numnodes, dof;
  fModel.CoordinateDimensions (numnodes, dof);
  fNodeMap.Allocate (numnodes);
  fModel.AllNodeMap (fNodeMap);
  fOutput->SetCoordinates (fModel.Coordinates(), &fNodeMap);
  cout << "\n Number of Nodes: " << numnodes << " DOF: " << dof << endl;
}

void TranslateIOManager::WriteNodeSets (void)
{
  int num = fModel.NumNodeSets ();
  if (num <= 0) return;
  ArrayT<StringT> names (num);
  fModel.NodeSetNames (names);

  int selection;
  cout << "\n Number of Node Sets: " << num << endl;
  cout << "\n1. Translate All\n";
  cout << "2. Translate Some\n";
  cout << "3. Translate None\n";
  cin >> selection;

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer ("no");
      if (selection == 2)
	{
	  cout << "    Translate Node Set " << names[i] << " (y/n) ? ";
	  cin >> answer;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y')
	fOutput->AddNodeSet (fModel.NodeSet (i), i+1);
    }
}
 
void TranslateIOManager::WriteElements (void)
{
  int num = fModel.NumElementGroups ();
  fOutputID.Allocate (num);
  ArrayT<StringT> names (num);
  fModel.ElementGroupNames (names);

  int selection;
  cout << "\n Number of Element Groups: " << num << endl;
  cout << "\n1. Translate All\n";
  cout << "2. Translate Some\n";
  cin >> selection;

  bool changing = false;
  StringT answer;
  for (int e=0; e < num; e++)
    {
      iArrayT block_ID (1);
      block_ID = e+1;

      // allow user to select element groups
      answer[0] = 'y';
      if (selection == 2)
	{
	  cout << "    Translate Element Group " << names[e] << " (y/n) ? ";
	  cin >> answer;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y')
	{
	  OutputSetT set (e, fModel.ElementGroupGeometry (e), block_ID, 
			  fModel.ElementGroup (e), fNodeLabels, fElementLabels, changing);
	  fOutputID[e] = fOutput->AddElementSet (set);
	}
      else
	fOutputID[e] = -1;
    }
}

void TranslateIOManager::WriteSideSets (void)
{
  int num = fModel.NumSideSets ();
  if (num <= 0) return;
  fGlobalSideSets.Allocate (num);
  ArrayT<StringT> names (num);
  fModel.SideSetNames (names);

  int selection;
  cout << "\n Number of Side Sets: " << num << endl;
  cout << "\n1. Translate All\n";
  cout << "2. Translate Some\n";
  cout << "3. Translate None\n";
  cin >> selection;

  if (selection == 3) return;
  for (int i=0; i < num; i++)
    {
      StringT answer ("yes");
      if (selection == 2)
	{
	  cout << "    Translate Node Set " << names[i] << " (y/n) ? ";
	  cin >> answer;
	}
      
      if (answer [0] == 'y' || answer[0] == 'Y')
	{
	  fGlobalSideSets[i] = fModel.SideSet (i);
	  int g = fModel.SideSetGroupIndex (i);
	  if (fModel.IsSideSetLocal (i))
	    fModel.SideSetLocalToGlobal (g, fGlobalSideSets[i], fGlobalSideSets[i]);
	  fOutput->AddSideSet (fGlobalSideSets [i], i+1, g);
	}
    }
}

void TranslateIOManager::VariableQuery (const ArrayT<StringT>& names, iArrayT& list)
{
  iAutoArrayT temp;
  for (int i=0; i < names.Length(); i++)
    {
      StringT answer;
      cout << " Extract variable " << names[i] << " (y/n) ? ";
      cin >> answer;

      if (answer[0] == 'y' || answer[0] == 'Y')
	temp.Append (i);
    }

  list.Allocate (temp.Length());
  list.CopyPart (0, temp, 0, temp.Length());  
}
