// file: MakeCSEFEManager.h

// MakeCSE version

// created: 11/10/99 SAW

#ifndef _MakeCSEFE_MANAGER_H_
#define _MakeCSEFE_MANAGER_H_

#include "ArrayT.h"
#include "MakeCSEIOManager.h"
#include "GlobalEdgeFinderT.h"
#include "NodeManagerPrimitive.h"
#include "MakeCSE.h"

using namespace Tahoe;

class MakeCSEFEManager
{
 public:
  MakeCSEFEManager (ostream& out, MakeCSEIOManager& theIO);
  ~MakeCSEFEManager (void);

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
  MakeCSE_ElementBaseT* ElementGroup(int groupnumber) const;

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
  ArrayT<MakeCSE_ElementBaseT*> fElementGroups;
  int fNumElementGroups;
  int fNumRegular;
  int fNumCSE;

  MakeCSE fCSEMakerBoss;
  GlobalEdgeFinderT fEdger;
};

inline NodeManagerPrimitive* MakeCSEFEManager::NodeManager(void) const { return fNodeBoss; }
inline int MakeCSEFEManager::NumCSEGroups (void) const { return fNumCSE; }
inline int MakeCSEFEManager::NumRegularGroups (void) const { return fNumRegular; }

#endif
