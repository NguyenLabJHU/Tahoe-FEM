
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
}

void ExtractIOManager::InitializeVariables (void)
{
  fNumNV = fModel.NumNodeVariables ();
  cout << "\n" << setw (10) << fNumNV << " Node Variables\n\n";
  fNodeLabels.Allocate (fNumNV);
  if (fNumNV > 0) fModel.NodeLabels (fNodeLabels);

  // query user as to which variables to translate
  VariableQuery (fNodeLabels, fNVUsed);
}

void ExtractIOManager::InitializeNodePoints (void)
{
  int selection;
  cout << "\n One file will be written per node point.\n";
  cout << "1. List of nodes\n";
  cout << "2. Node Set\n";
  cout << "3. Every nth node\n";
  cout << "\n How do you want to define your list of nodes: ";
  cin >> selection;

  int numnodes, numdims;
  fModel.CoordinateDimensions (numnodes, numdims);
  fNodeMap.Allocate (numnodes);
  fModel.AllNodeMap (fNodeMap);
  switch (selection)
    {
    case 1:
      {
	cout << "\n Enter the number of nodes: ";
	cin >> fNumNP;
	fNodePoints.Allocate (fNumNP);
	fNodePointIndex.Allocate (fNumNP);
	for (int n=0; n < fNumNP; n++)
	  {
	    cout << " Enter node " << n+1 << ": ";
	    cin >> fNodePoints[n];

	    // translate node numbers to index
	    int index;
	    fNodeMap.HasValue (fNodePoints[n], index);
	    if (index < 0 || index >= numnodes) 
	      {
		cout << " ExtractIOManager::InitializeNodePoints\n";
		cout << " Node " << fNodePoints[n] << " was not found.\n";
		throw eOutOfRange;
	      }
	    fNodePointIndex [n] = index;
	  }
	break;
      }
    case 2:
      {
	int num = fModel.NumNodeSets ();
	ArrayT<StringT> nodesetnames (num);
	fModel.NodeSetNames (nodesetnames);
	cout << "\n";
	for (int h=0; h < num; h++)
	  cout << "    " << h+1 << ". " << nodesetnames[h] << "\n";
	cout << "\n Enter the number of the node set: ";
	int ni;
	cin >> ni;
	ni--;
	fNumNP = fModel.NodeSetLength (ni);
	fNodePoints.Allocate (fNumNP);
	fNodePointIndex.Allocate (fNumNP);
	fNodePointIndex = fModel.NodeSet (ni);
	for (int n=0; n < fNumNP; n++)
	  fNodePoints[n] = fNodeMap [fNodePointIndex[n]];
	break;
      }
    case 3:
      {
	int freq;
	cout << "\n Number of Nodes: " << numnodes << "\n";
	cout << "   Enter n: ";
	cin >> freq;
	fNumNP = numnodes/freq;
	fNodePoints.Allocate (fNumNP);
	fNodePointIndex.Allocate (fNumNP);
	for (int n=0; n < fNumNP; n++)
	  {
	    fNodePoints[n] = fNodeMap[n*freq];
	    fNodePointIndex[n] = n*freq;
	  }
	break;
      }
    default:
      throw eGeneralFail;
    }
}

void ExtractIOManager::TranslateVariables (void)
{
  // open files, one per node
  StringT temp, ext;
  temp.Append (fNodePoints[fNumNP-1]);
  int digits = temp.Length()-1;
  switch (fOutputFormat)
    {
    case IOBaseT::kTecPlot: ext = "dat"; break;
    case IOBaseT::kTahoe: ext = "txt"; break;
    }
  PrepFiles (ext, digits);

  // keep user informed
  int check = 1;
  if (fNumTS > 100) check = 10;
  else if (fNumTS > 1000) check = 100;
  else if (fNumTS > 10000) check = 1000;

  // read data
  int numnodes, numdims;
  int numused = fNVUsed.Length();
  fModel.CoordinateDimensions (numnodes, numdims);
  dArray2DT nv (numnodes, fNumNV);
  for (int t=0; t < fNumTS; t++)
    {
      if ((t+1)%check == 0 || t == fNumTS-1)
	cout << "Time Step " << fTimeIncs[t]+1 << ": " << fTimeSteps[t] << endl;
      fModel.AllNodeVariables (fTimeIncs[t], nv);

      // write files
      switch (fOutputFormat)
	{
	case IOBaseT::kTecPlot: // point format
	case IOBaseT::kTahoe:
	  {
	    for (int n=0; n < fNumNP; n++)
	      {
		ofstreamT outfile;
		OpenFile (outfile, n, digits, ext, true);
		outfile << fTimeSteps[t] << " ";
		for (int g=0; g < numused; g++)
		  outfile << nv (fNodePointIndex[n], fNVUsed[g]) << " ";
		outfile << "\n";
	      }
	    break;
	  }
	default:
	  throw eGeneralFail;
	}
    }
}

void ExtractIOManager::PrepFiles (StringT& ext, int digits) const
{
  for (int i=0; i < fNumNP; i++)
    {
      ofstreamT outfile;
      OpenFile (outfile, i, digits, ext, false);
      switch (fOutputFormat)
	{
	case IOBaseT::kTecPlot:
	  {
	    iArrayT ijk (2);
	    ijk[0] = fNumTS;
	    ijk[1] = 1;

	    ArrayT<StringT> labels (fNVUsed.Length() + 1);
	    labels[0] = "Time";
	    for (int il=0; il < fNVUsed.Length(); il++)
	      labels[il+1] = fNodeLabels [fNVUsed[il]];

	    TecPlotT tec (fMessage, true);
	    tec.WriteHeader (outfile, outfile.filename(), labels);
	    tec.WriteIJKZone (outfile, outfile.filename(), ijk);
	    break;
	  }
	}
    }
}

void ExtractIOManager::OpenFile (ofstreamT& o, int index, int digits, StringT& ext, bool append) const
{
  StringT filename (fOutputName);
  filename.Append ("_", fNodePoints[index], digits);
  filename.Append (".", ext);
  if (!append)
    {
      remove (filename); // remove any pre-existing file
      o.open (filename);
    }
  else
    o.open_append (filename);
  if (!o.is_open())
    {
      fMessage << "\nExtractIOManager::OpenFile cannot open file: "
	       << filename << "\n\n";
      throw eGeneralFail;
    }
}
