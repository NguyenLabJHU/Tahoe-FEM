// file: FEManager.cpp

// MakeCSE version

// created: 11/10/99 SAW

#include "FEManager.h"
#include "Quad2Tri.h"
#include "CSEBaseT.h"

using namespace Tahoe;

FEManager::FEManager (ostream& out, MakeCSEIOManager& theIO) :
  fMainOut (out),
  fEdger (out),
  fCSEMakerBoss (out, fEdger)
{
  theIO.InputData (fPrintInput, MakeCSEIOManager::kExecution);

  SetNodeManager (theIO);
  SetElementGroups (theIO);

  fEdger.Initialize (*this, fNodeBoss->NumNodes());
  fCSEMakerBoss.Initialize (theIO, *this, fPrintInput);

  theIO.InputData (fRenumberOption, MakeCSEIOManager::kRenumber);
  if (fRenumberOption < kNoRenumber || fRenumberOption > kRenumberAll) 
    fRenumberOption = kNoRenumber;
}

FEManager::~FEManager (void)
{
  delete fNodeBoss;
  for (int i=0; i < fElementGroups.Length(); i++)
    delete fElementGroups[i];
}

void FEManager::CreateCSE (void)
{
  fCSEMakerBoss.Create ();

  if (fRenumberOption != kNoRenumber)
    {
      iArrayT map;
      fNodeBoss->Renumber(fRenumberOption, map);

      cout << "\n Renumbering Connectivity " << endl;
      for (int e=0; e < fElementGroups.Length(); e++)
	fElementGroups[e]->Renumber (map);
    }
}

ElementBaseT* FEManager::ElementGroup(int groupnumber) const
{
  // check range 
	if (groupnumber > -1 && groupnumber < fNumElementGroups)
		return fElementGroups[groupnumber];
	else
		return NULL;
}

void FEManager::NodesUsed (int groupID, iArrayT& nodes) const
{
  int g = -1;
  for (int i=0; i < fElementGroups.Length(); i++)
    if (fElementGroups[i]->GroupNumber() == groupID) g = i;

  nodes.Allocate (0);
  if (g > -1) fElementGroups[g]->NodesUsed (nodes);
}

void FEManager::SetIO (MakeCSEIOManager& theIO)
{
  fNodeBoss->RegisterOutput (theIO);

  for (int i=0; i < fElementGroups.Length(); i++)
    fElementGroups[i]->RegisterOutput (theIO);
}

void FEManager::WriteOutput (MakeCSEIOManager& theIO, IOBaseT::OutputModeT mode) const
{
  for (int i = 0; i < fElementGroups.Length(); i++)
    fElementGroups[i]->WriteOutput (theIO, mode);
}

//************** PRIVATE *******************

void FEManager::SetNodeManager (MakeCSEIOManager& theIO)
{
  fMainOut << "\n N o d a l   D a t a :\n\n";
  fNodeBoss = new NodeManagerPrimitive (fMainOut, fPrintInput, *this);
  fNodeBoss->Initialize (theIO);
}

void FEManager::SetElementGroups (MakeCSEIOManager& theIO)
{
  // determine number of regular element groups 
  fNumRegular = theIO.NumElementGroups ();
  const ArrayT<StringT> sids = theIO.ElementGroupIDs ();
  iArrayT ids (sids.Length());
  for (int g=0; g < sids.Length(); g++)
    ids[g] = atoi (sids[g]);

  // determine number of ppossible CSE element groups
  iArrayT facets, zones, boundaries;
  theIO.InputData (facets, MakeCSEIOManager::kFacet);
  theIO.InputData (zones, MakeCSEIOManager::kZone);
  theIO.InputData (boundaries, MakeCSEIOManager::kBoundary);
  iAutoArrayT CSEids;
  for (int f=0; f < facets.Length(); f += 3)
    CSEids.AppendUnique (facets [f + 2]);
  for (int z=0; z < zones.Length(); z += 2)
    CSEids.AppendUnique (zones [z + 1]);
  for (int b=0; b < boundaries.Length(); b += 4)
    CSEids.AppendUnique (boundaries [b + 3]);
  fNumCSE = CSEids.Length();

  // print data
  fMainOut << "\n E l e m e n t   D a t a :\n\n";
  fMainOut << " Number of regular element groups. . . . . . . . = "
	   << fNumRegular << '\n';
  ids.WriteWrapped (fMainOut, 6);
  fMainOut << " Number of CSE element group . . . . . . . . . . = " 
	   << fNumCSE << '\n' << CSEids;

  // see if element groups will be split
  iArrayT split;
  theIO.InputData (split, MakeCSEIOManager::kElementSplit);
  iArray2DT splitData;
  splitData.Set (split.Length()/2, 2, split.Pointer());

  // set up and initialize regular element groups
  fElementGroups.Allocate (fNumRegular + fNumCSE);
  fNumElementGroups = fNumRegular + fNumCSE;
  for (int e=0; e < fNumRegular; e++)
    {
      // see if element group is to be split
      int dex;
      if (splitData.ColumnHasValue (0, ids[e], dex))
	{
	  switch (splitData(dex, 1))
	    {
	    case Quad2Tri::kXMethod:
	    case Quad2Tri::kSlashMethod:
	    case Quad2Tri::kBackSlashMethod:
	    case Quad2Tri::kStarMethod:
	      fElementGroups[e] = new Quad2Tri (fMainOut, *fNodeBoss, splitData (dex, 1), ids[e]);
	      break;
	    default:
	      cout << "FEManager::SetElementGroups: unable to split element block: " << ids[e] << ", unknown method: " << splitData (dex, 1) << endl;
	      throw eGeneralFail;
	    }
	}

      // if not split
      else
	fElementGroups[e] = new ElementBaseT (fMainOut, ids[e]);

      fElementGroups[e]->Initialize (theIO);
    }

  // set up CSE element groups, finish initialize during MakeCSE
  for (int c=0; c < fNumCSE; c++)
    fElementGroups[c + fNumRegular] = new CSEBaseT (fMainOut, CSEids[c]);
}

