
#include "ExtractIOManager.h"
#include "TecPlotT.h"

ExtractIOManager::ExtractIOManager (ostream& out) :
  TranslateIOManager (out)
{
}

void ExtractIOManager::Translate (const StringT& program, const StringT& version, const StringT& title)
{
  fModel.Initialize ();
  SetOutput (program, version, title);

  InitializeVariables ();
  if (fNumNV < 1)
    {
      fMessage << "\n No nodal variables found.";
      return;
    }

  InitializeTime();
  if (fNumTS < 1)
    {
      fMessage << "\n No time steps found.";
      return;
    }

  InitializeNodePoints();
  if (fNumNP < 1)
    {
      fMessage << "\n No node points found.";
      return;
    }

  TranslateVariables ();
}

/**************** PRIVATE **********************/

void ExtractIOManager::SetOutput (const StringT& program, const StringT& version, const StringT& title)
{
  cout << "\n    eq.  " << IOBaseT::kTahoe   << ". Text\n";
  cout << "    eq.  " << IOBaseT::kTecPlot << ". TecPlot 7.5\n";
  cout << "\n Enter the Output Format: ";
  cin >> fOutputFormat;
  cout << "\n Enter the root of the output files: ";
  cin >> fOutputName;

  StringT ext;
  switch (fOutputFormat)
    {
    case IOBaseT::kTecPlot:
      ext = ".dat";
      break;
    case IOBaseT::kTahoe:
      ext = ".txt";
      break;
    }  
  StringT filename (fOutputName);
  filename.Append (ext);
  fOutFile.open (filename);
}

void ExtractIOManager::InitializeVariables (void)
{
  fNumNV = fModel.NumNodeVariables ();
  cout << "\n" << setw (10) << fNumNV << " Node Variables\n";
  fNodeLabels.Allocate (fNumNV);
  if (fNumNV > 0) fModel.NodeLabels (fNodeLabels);
  // future: query user as to which variables to translate
}

void ExtractIOManager::InitializeNodePoints (void)
{
  fNumNP = 1;
  fNodePoints.Allocate (fNumNP);
  fNodePointIndex.Allocate (fNumNP);
  cout << "\n Enter the node number to extract data at: ";
  cin >> fNodePoints[0];

  int numnodes, numdims;
  fModel.CoordinateDimensions (numnodes, numdims);
  fNodeMap.Allocate (numnodes);
  fModel.AllNodeMap (fNodeMap);

  int index;
  fNodeMap.HasValue (fNodePoints[0], index);
  if (index < 0 || index >= numnodes) throw eOutOfRange;
  fNodePointIndex [0] = index;
}

void ExtractIOManager::TranslateVariables (void)
{
  dArray2DT values (fNumTS, fNumNV);

  int numnodes, numdims;
  fModel.CoordinateDimensions (numnodes, numdims);
  dArray2DT nv (numnodes, fNumNV);

  int check = 1;
  if (fNumTS > 100) check = 10;
  if (fNumTS > 1000) check = 100;
  if (fNumTS > 10000) check = 1000;
  for (int t=0, offset= 0; t < fNumTS; t++, offset += fNumNV)
    {
      if ((t+1)%check == 0 || t == fNumTS-1)
	cout << "Time Step " << fTimeIncs[t]+1 << ": " << fTimeSteps[t] << endl;
      fModel.AllNodeVariables (fTimeIncs[t], nv);
      values.CopyPart (offset, nv, fNodePointIndex[0]*fNumNV, fNumNV);
    }

  switch (fOutputFormat)
    {
    case IOBaseT::kTecPlot:
      {
	StringT nodenumber;
	nodenumber.Append (fNodePoints[0]);
	iArrayT ijk (2);
	ijk[0] = fNumTS;
	ijk[1] = 1;
	ArrayT<StringT> labels (fNumNV + 1);
	labels[0] = "Time";
	for (int i=0; i < fNumNV; i++)
	  labels[i+1] = fNodeLabels[i];

	TecPlotT tec (fMessage, false);
	tec.WriteHeader (fOutFile, fOutputName, labels);
	tec.WriteIJKZone (fOutFile, nodenumber, ijk);
	tec.WriteData (fOutFile, fTimeSteps, fNumTS, 1);
	tec.WriteData (fOutFile, values, fNumTS, fNumNV);
	break;
      }
    case IOBaseT::kTahoe:
      {
	for (int i=0; i < fNumTS; i++)
	  {
	    fOutFile << fTimeSteps[i] << " ";
	    values.PrintRow (i, fOutFile);
	  }
      }
    default:
      return;
    }
  fOutFile.close();
}
