/* $Id: TranslateIOManager.cpp,v 1.2 2001-09-07 13:25:29 sawimme Exp $  */

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
  SetOutput (program, version, program);

  InitializeVariables ();
  if (fNumNV < 1 && fNumEV < 1) 
    {
      fMessage << "\n No variables found, writing geometry only.";
      WriteGeometry ();
      fOutput->WriteGeometry ();
      return;
    }
  
  StringT answer;
  bool xy = false;
  // future: extract xy data
  //cout << "\n Do you want variables extracted as xy data only (y/n)? \n";
  //cin >> answer;
  //if (answer[0] == 'y' || answer[0] == 'Y') 
  //xy = true;

  if (!xy)
    WriteGeometry ();
  //else
  // open Tecplot or text or Matlab ???

  InitializeTime ();
  if (fNumTS < 1)
    {
      fMessage << "\n No time steps found, writing geometry only.";
      WriteGeometry ();
      fOutput->WriteGeometry ();
      return;      
    }
  
  TranslateVariables (xy);
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

  fMessage << "\n" << setw (10) << fNumNV << " Node Variables\n";
  fMessage << setw (10) << fNumEV << " Element Variables\n";
  fMessage << setw (10) << fNumQV << " Quadrature Varaibles\n";

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
    fModel.TimeSteps (fTimeSteps);

  // future: query user as to which time steps to translate
}

void TranslateIOManager::TranslateVariables (bool xy)
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

  for (int t=0; t < fNumTS; t++)
    {
      // future: see if user wanted to translate this time step

      fMessage << "Time Step " << t+1 << ": " << fTimeSteps[t] << endl;
      for (int g=0; g < fOutputID.Length(); g++)
	{
	  // read all variables for this step
	  dArray2DT nvalues (0,0), evalues (0,0);
	  if (fNumNV > 0) 
	    fModel.NodeVariables (t, groupnames[g], nvalues);
	  if (fNumEV > 0)
	    fModel.ElementVariables (t, groupnames[g], evalues);

	  // future: accout for which variables user wanted translated

	  if (!xy)
	    fOutput->WriteOutput (fTimeSteps[t], fOutputID[g], nvalues, evalues);
	  //else
	  // future: extract xy data
	}
    }
}

void TranslateIOManager::WriteGeometry (void)
{
  fMessage << "\n";
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
  fMessage << setw (10) << numnodes << " Nodes\n";
  fMessage << setw (10) << dof << " DOF\n";
}

void TranslateIOManager::WriteNodeSets (void)
{
  int num = fModel.NumNodeSets ();
  for (int i=0; i < num; i++)
    fOutput->AddNodeSet (fModel.NodeSet (i), i+1);
  fMessage << setw (10) << num << " Node Sets\n";
}
 
void TranslateIOManager::WriteElements (void)
{
  int num = fModel.NumElementGroups ();
  fOutputID.Allocate (num);
  bool changing = false;
  for (int e=0; e < num; e++)
    {
      iArrayT block_ID (1);
      block_ID = e+1;
      //cout << e << " " << ElementGroupGeometry(e) << endl;
      
      OutputSetT set (e, fModel.ElementGroupGeometry (e), block_ID, 
		      fModel.ElementGroup (e),
		      fNodeLabels, fElementLabels, changing);
      fOutputID[e] = fOutput->AddElementSet (set);
    }
  fMessage << setw (10) << num << " Element Groups\n";
}

void TranslateIOManager::WriteSideSets (void)
{
  int num = fModel.NumSideSets ();
  fGlobalSideSets.Allocate (num);
  for (int i=0; i < num; i++)
    {
      fGlobalSideSets[i] = fModel.SideSet (i);
      int g = fModel.SideSetGroupIndex (i);
      if (fModel.IsSideSetLocal (i))
	fModel.SideSetLocalToGlobal (g, fGlobalSideSets[i], fGlobalSideSets[i]);
      fOutput->AddSideSet (fGlobalSideSets [i], i+1, g);
    }
  fMessage << setw (10) << num << " Side Sets\n";
}
