// $Id: MakeCSE_FEManager.h,v 1.4 2003-09-05 23:11:47 paklein Exp $
// created: 11/10/99 SAW
#ifndef _MakeCSE_FE_MANAGER_H_
#define _MakeCSE_FE_MANAGER_H_

#include "ArrayT.h"
#include "MakeCSE_IOManager.h"
#include "GlobalEdgeFinderT.h"
#include "NodeManagerPrimitive.h"
#include "MakeCSE.h"
#include "CSEConstants.h"
#include "OutputBaseT.h"
#include "ModelManagerT.h"

namespace Tahoe {

class MakeCSE_FEManager
{
 public:
  MakeCSE_FEManager (ostream& out);
  ~MakeCSE_FEManager (void);

  void InitializeInput (ifstreamT& in, bool interactive);
  void InitializeOutput (const StringT& title, const StringT& program_name, const StringT& version);
  void CreateCSE (void);

  /* accessors */
  int NumRegularGroups (void) const;
  int NumCSEGroups (void) const;
  NodeManagerPrimitive* NodeManager(void) const;
  MakeCSE_ElementBaseT* ElementGroup(int groupnumber) const;

  void NodesUsed (const StringT& groupID, iArrayT& nodes) const;

  void WriteOutput (void);

 private:
  void SetNodeManager (void);
  void SetElementGroups (void);
  void ReadCSEIDs (sArrayT& cseids);

 private:
  ostream& fMainOut;
  bool fPrintInput;

  MakeCSE_IOManager*   fParameters;
  ModelManagerT        fModel;
  OutputBaseT*         fOutput;

  NodeManagerPrimitive* fNodeBoss;
  ArrayT<MakeCSE_ElementBaseT*> fElementGroups;
  int fNumElementGroups;
  int fNumRegular;
  int fNumCSE;

  MakeCSE fCSEMakerBoss;
  GlobalEdgeFinderT fEdger;
};

inline NodeManagerPrimitive* MakeCSE_FEManager::NodeManager(void) const { return fNodeBoss; }
inline int MakeCSE_FEManager::NumCSEGroups (void) const { return fNumCSE; }
inline int MakeCSE_FEManager::NumRegularGroups (void) const { return fNumRegular; }
}
#endif
