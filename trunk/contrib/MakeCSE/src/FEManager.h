// file: FEManager.h

// MakeCSE version

// created: 11/10/99 SAW

#ifndef _FE_MANAGER_H_
#define _FE_MANAGER_H_

#include "ArrayT.h"
#include "MakeCSEIOManager.h"
#include "GlobalEdgeFinderT.h"
#include "NodeManagerPrimitive.h"
#include "MakeCSE.h"

using namespace Tahoe;

class FEManager
{
 public:
  FEManager (ostream& out, MakeCSEIOManager& theIO);
  ~FEManager (void);

  void CreateCSE (void);

  enum settings {kNotSet = GlobalEdgeFinderT::kNoNeighbor,
		 kExteriorFacet = GlobalEdgeFinderT::kExteriorFacet};

  enum renumberoptions { kNoRenumber = 0,
			 kRenumberAdded = 1,
			 kRenumberAll = 2 };

  /* accessors */
  int NumRegularGroups (void) const;
  int NumCSEGroups (void) const;
  NodeManagerPrimitive* NodeManager(void) const;
  ElementBaseT* ElementGroup(int groupnumber) const;

  void NodesUsed (int groupID, iArrayT& nodes) const;

  void SetIO (MakeCSEIOManager& theIO);
  void WriteOutput (MakeCSEIOManager& theIO, IOBaseT::OutputModeT mode) const;

 private:
  void SetNodeManager (MakeCSEIOManager& theIO);
  void SetElementGroups (MakeCSEIOManager& theIO);

 private:
  ostream& fMainOut;
  int fPrintInput;
  int fRenumberOption;

  NodeManagerPrimitive* fNodeBoss;
  ArrayT<ElementBaseT*> fElementGroups;
  int fNumElementGroups;
  int fNumRegular;
  int fNumCSE;

  MakeCSE fCSEMakerBoss;
  GlobalEdgeFinderT fEdger;
};

inline NodeManagerPrimitive* FEManager::NodeManager(void) const { return fNodeBoss; }
inline int FEManager::NumCSEGroups (void) const { return fNumCSE; }
inline int FEManager::NumRegularGroups (void) const { return fNumRegular; }

#endif
