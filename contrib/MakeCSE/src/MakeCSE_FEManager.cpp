// file: FEManager.cpp

// MakeCSE version

// created: 11/10/99 SAW

#include "MakeCSE_FEManager.h"
#include "Quad2Tri.h"
#include "MakeCSE_CSEBaseT.h"
#include "InteractiveIOManagerT.h"
#include "ParameterFileManagerT.h"
#include "ifstreamT.h"
#include "sArrayT.h"

using namespace Tahoe;

MakeCSE_FEManager::MakeCSE_FEManager (ostream& out) :
  fMainOut (out),
  fEdger (out),
  fCSEMakerBoss (out, fEdger),
  fModel (out),
  fParameters (NULL)
{
}

MakeCSE_FEManager::~MakeCSE_FEManager (void)
{
  delete fNodeBoss;
  for (int i=0; i < fElementGroups.Length(); i++)
    delete fElementGroups[i];
}

void MakeCSE_FEManager::InitializeInput (ifstreamT& in, bool interactive)
{
  IOBaseT::FileTypeT format;
  StringT name;
  if (interactive)
    {
      fMainOut << " No input file, interactive session.\n\n";
      fParameters = new InteractiveIOManagerT ();
      fParameters->Initialize ();
      fModel.Initialize ();
      fModel.Format(format, name);
      fParameters->InputFormat (format, name); // echo the format
    }
  else 
    {
      fMainOut << " Reading from Input file . . . . . . . . . . . . = " 
	  << in.filename() << '\n';
      fParameters = new ParameterFileManagerT (in.filename());
      fParameters->Initialize ();
      fParameters->InputFormat (format, name);
      IOBaseT temp (fMainOut);
      format = temp.int_to_FileTypeT (format);
      fModel.Initialize (format, name, true);
    }

  fModel.EchoData (fMainOut);
  fPrintInput = fParameters->Verbose();

  SetNodeManager ();
  SetElementGroups ();

  fEdger.Initialize (*this, fNodeBoss->NumNodes());
  fCSEMakerBoss.Initialize (fModel, *fParameters, *this, fPrintInput);
}

void MakeCSE_FEManager::CreateCSE (void)
{
  fCSEMakerBoss.Create ();

  CSEConstants::RenumberMethodT meth = fParameters->RenumberMethod();
  if (meth != CSEConstants::kNoRenumber)
    {
      iArrayT map;
      fNodeBoss->Renumber(meth, map);

      cout << "\n Renumbering Connectivity " << endl;
      for (int e=0; e < fElementGroups.Length(); e++)
	fElementGroups[e]->Renumber (map);
    }
}

MakeCSE_ElementBaseT* MakeCSE_FEManager::ElementGroup(int groupnumber) const
{
  // check range 
	if (groupnumber > -1 && groupnumber < fNumElementGroups)
		return fElementGroups[groupnumber];
	else
		return NULL;
}

void MakeCSE_FEManager::NodesUsed (const StringT& groupID, iArrayT& nodes) const
{
  int g = -1;
  for (int i=0; i < fElementGroups.Length(); i++)
    if (fElementGroups[i]->GroupNumber() == groupID) g = i;

  nodes.Allocate (0);
  if (g > -1) fElementGroups[g]->NodesUsed (nodes);
}

void MakeCSE_FEManager::WriteOutput (void) const
{
}

//************** PRIVATE *******************

void MakeCSE_FEManager::SetNodeManager (void)
{
  fMainOut << "\n N o d a l   D a t a :\n\n";
  fNodeBoss = new NodeManagerPrimitive (fMainOut, fPrintInput, *this);
  fNodeBoss->Initialize (fModel, *fParameters);
}

void MakeCSE_FEManager::SetElementGroups (void)
{
  // determine number of regular element groups 
  fNumRegular = fModel.NumElementGroups ();
  const sArrayT ids = fModel.ElementGroupIDs ();

  // determine number of ppossible CSE element groups
  sArrayT facets, zones, boundaries;
  fParameters->Facets (facets);
  fParameters->Zones (zones);
  fParameters->Boundaries (boundaries);

  fNumCSE = facets.Length()/3 + zones.Length()/2 + boundaries.Length()/3;
  sArrayT CSEids (fNumCSE);
  int count = 0;
  for (int f=0; f < facets.Length(); f += 3)
    CSEids[count++] = facets [f + 2];
  for (int z=0; z < zones.Length(); z += 2)
    CSEids[count++] = zones [z + 1];
  for (int b=0; b < boundaries.Length(); b += 4)
    CSEids[count++] = boundaries [b + 3];

  // print data
  fMainOut << "\n E l e m e n t   D a t a :\n\n";
  fMainOut << " Number of regular element groups. . . . . . . . = "
	   << fNumRegular << '\n';
  ids.WriteWrapped (fMainOut, 4);
  fMainOut << '\n';
  fMainOut << " Number of CSE element group . . . . . . . . . . = " 
	   << fNumCSE << '\n';
  CSEids.WriteWrapped (fMainOut, 4);
  fMainOut << '\n';

  // see if element groups will be split
  sArrayT split;
  ArrayT<CSEConstants::SplitMethodT> methods;
  fParameters->SplitBlocks (split, methods);

  // set up and initialize regular element groups
  fElementGroups.Allocate (fNumRegular + fNumCSE);
  fNumElementGroups = fNumRegular + fNumCSE;
  for (int e=0; e < fNumRegular; e++)
    {
      // see if element group is to be split
      int dex = -1;
      if (split.HasValue (ids[e], dex))
	{
	  switch (methods[dex])
	    {
	    case CSEConstants::kXMethod:
	    case CSEConstants::kSlashMethod:
	    case CSEConstants::kBackSlashMethod:
	    case CSEConstants::kStarMethod:
	      fElementGroups[e] = new Quad2Tri (fMainOut, *fNodeBoss, methods[dex], ids[e]);
	      break;
	    default:
	      cout << "MakeCSE_FEManager::SetElementGroups: unable to split element block: " << ids[e] << ", unknown method: " << methods[dex] << endl;
	      throw eGeneralFail;
	    }
	}

      // if not split
      else
	fElementGroups[e] = new MakeCSE_ElementBaseT (fMainOut, ids[e]);

      fElementGroups[e]->Initialize (fModel, *fParameters);
    }

  // set up CSE element groups, finish initialize during MakeCSE
  for (int c=0; c < fNumCSE; c++)
    fElementGroups[c + fNumRegular] = new MakeCSE_CSEBaseT (fMainOut, CSEids[c]);
}

