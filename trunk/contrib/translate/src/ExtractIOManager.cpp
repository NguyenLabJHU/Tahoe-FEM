/* $Id: ExtractIOManager.cpp,v 1.17 2004-12-13 23:47:08 paklein Exp $ */
#include "ExtractIOManager.h"

#include "ExceptionT.h"
#include "TecPlotT.h"
#include <stdio.h>

using namespace Tahoe;

ExtractIOManager::ExtractIOManager (ostream& out, istream& in, bool write) :
  TranslateIOManager (out, in, write)
{
}

void ExtractIOManager::Translate (const StringT& program, const StringT& version, const StringT& title)
{
  // set input
  SetInput ();
  
  cout << "\n" << setw (10) << fModel.NumNodeVariables () << " Node Variables\n";
  cout << setw (10) << fModel.NumElementVariables () << " Element Variables\n";
  cout << setw (10) << fModel.NumQuadratureVariables () << " Quadrature Variables\n";

  // set output
  SetOutput (program, version, title);

  // set time steps
  InitializeTime();
  if (fNumTS < 1)
    {
      fMessage << "\n No time steps found.";
      return;
    }
  // keep user informed
  if (fNumTS < 10) fCheck = fNumTS;
  else if (fNumTS < 100) fCheck = int(fNumTS * 0.1);
  else if (fNumTS < 1000) fCheck = int(fNumTS * 0.05);
  else fCheck = 100;

  // set variables
  Initialize ();

  // go
  TranslateVariables ();
}

/**************** PROTECTED **********************/

void ExtractIOManager::SetOutput (const StringT& program, const StringT& version, const StringT& title)
{
  fOutputFormat = -1;
  while (fOutputFormat != IOBaseT::kTahoe &&
	 fOutputFormat != IOBaseT::kTecPlot)
    {
      if (fWrite)
	{
	  cout << "\n    eq.  " << IOBaseT::kTahoe   << ". Text\n";
	  cout << "    eq.  " << IOBaseT::kTecPlot << ". TecPlot 7.5\n";
	  cout << "\n Enter the Output Format: ";
	}
      fIn >> fOutputFormat;
      if (fEcho) fEchoOut << fOutputFormat << "\n";
    }
  if (fWrite)
    cout << "\n Enter the root of the output files: ";
  fIn >> fOutputName;
  if (fEcho) fEchoOut << fOutputName << "\n";

  switch (fOutputFormat)
    {
    case IOBaseT::kTecPlot: fOutfileExtension = "dat"; break;
    case IOBaseT::kTahoe: fOutfileExtension = "txt"; break;
    }
  cout << "\n Output Format: " << fOutputFormat << " File: " << fOutputName << endl;
}

void ExtractIOManager::PrepFiles (iArrayT& varsused, ArrayT<StringT>& labels)
{
  cout << "\n";

  // open files, one per item = node or element or int point
  fNumDigits = fItemNames[0].StringLength();
  for (int d=1; d < fItemNames.Length(); d++)
    if (fItemNames[d].StringLength() > fNumDigits) 
      fNumDigits = fItemNames[d].StringLength();

  int numused = varsused.Length() + 1; // account for time variable
  if (fCoords > 0) // account for coordinates as a variable
    numused += fCoords;

  for (int i=0; i < fNumItems; i++)
    {
      ofstreamT outfile;
      OpenFile (outfile, fItemNames[i], false);
      switch (fOutputFormat)
	{
	case IOBaseT::kTecPlot:
	  {
	    iArrayT ijk (2);
	    ijk[0] = fNumTS;
	    ijk[1] = 1;

	    ArrayT<StringT> printlabels (numused);
	    printlabels[0] = "Time";
	    for (int ic=0; ic < fCoords; ic++)
	      {
		printlabels[ic + 1] = "X";
		printlabels[ic + 1].Append (ic+1);
	      }
	    for (int il=0; il < varsused.Length(); il++)
	      printlabels[il + 1 + fCoords] = labels [varsused[il]];

	    TecPlotT tec (fMessage, true);
	    tec.WriteHeader (outfile, outfile.filename(), printlabels);
	    tec.WriteIJKZone (outfile, fItemNames[i], ijk);
	    break;
	  }
	}
    }
}

void ExtractIOManager::WriteVarData (iArrayT& varsused, int ts) const
{
  if ((ts+1)%fCheck == 0 || ts == fNumTS-1)
    {
      /*time_t t;
	time (&t);
	cout << " STATUS:  Time Step " << fTimeIncs[ts]+1 << ": " << fTimeSteps[ts] 
	<< " Wall Time: " << localtime (&t) << endl; */
      cout << " STATUS:  Time Step " << fTimeIncs[ts]+1 << ": " << fTimeSteps[ts] << endl;
    }

  switch (fOutputFormat)
    {
    case IOBaseT::kTecPlot: // point format
    case IOBaseT::kTahoe:
      {
	for (int n=0; n < fNumItems; n++)
	  {
	    ofstreamT outfile;
	    OpenFile (outfile, fItemNames[n], true);

	    outfile << fTimeSteps[ts] << " ";

	    for (int ic=0; ic < fCoords; ic++)
	      outfile << fCoordinates (fItemIndex[n], ic) << " ";

	    for (int g=0; g < varsused.Length(); g++)
	      outfile << fVarData (fItemIndex[n], varsused[g]) << " ";
	    outfile << "\n";
	  }
	break;
      }
    default:
      ExceptionT::BadInputValue ("ExtractIOManager::WriteVarData","Unacceptable output format %i", fOutputFormat);
    }
}

/**************** PRIVATE **********************/

void ExtractIOManager::OpenFile (ofstreamT& o, const StringT& name, bool append) const
{
  StringT filename (fOutputName);
  filename.Append ("_");
  for (int i=name.StringLength(); i < fNumDigits; i++)
    filename.Append ("0");
  filename.Append (name);
  filename.Append (".", fOutfileExtension);
  if (!append)
    {
      remove (filename); // remove any pre-existing file
      o.open (filename);
    }
  else
    o.open_append (filename);
  if (!o.is_open())
    ExceptionT::GeneralFail ("ExtractIOManager::OpenFile","Unable to open %s", filename.Pointer());
}

