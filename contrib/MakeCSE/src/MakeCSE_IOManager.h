#ifndef _MAKECSEIOMANAGER_T_H_
#define _MAKECSEIOMANAGER_T_H_

#include "ios_fwd_decl.h"
#include "OutputBaseT.h"
#include "ofstreamT.h"
#include "CSEConstants.h"
#include "sArrayT.h"

namespace Tahoe {

/* forward declarations */
class ModelManagerT;

/** this class is in charge of reading parameter data from file or interactively */
class MakeCSE_IOManager
{
 public:
  MakeCSE_IOManager (void);


  virtual void Initialize (void) = 0;

  virtual void InputFormat (IOBaseT::FileTypeT &format, StringT& name) = 0;
  virtual void OutputFormat (IOBaseT::FileTypeT &format, StringT& name) = 0;
  virtual bool Verbose (void) = 0;
  virtual void Facets (sArrayT& names) = 0;
  virtual void Zones (sArrayT& names) = 0;
  virtual void Boundaries (sArrayT& names) = 0;
  virtual CSEConstants::ZoneEdgeT ZoneMethod (void) = 0;
  virtual void ZoneEdgeNodeSets (sArrayT& names) = 0;
  virtual void Contact (sArrayT& names) = 0;
  virtual void SingleNodes (sArrayT& names) = 0;
  virtual void BlockToNode (sArrayT& names) = 0;
  virtual void NodeSetsMapped (sArrayT& names, ArrayT<CSEConstants::NodeMapMethodT>& meths) = 0;
  virtual void SideSetsMapped (sArrayT& names) = 0;
  virtual CSEConstants::RenumberMethodT RenumberMethod (void) = 0;
  virtual void SplitBlocks (sArrayT& names, ArrayT<CSEConstants::SplitMethodT>& meths) = 0;

 protected:
  CSEConstants::CSEMethodT int2CSEMethodT (int i);
  CSEConstants::ZoneEdgeT int2ZoneEdgeT (int i);
  CSEConstants::NodeMapMethodT int2NodeMapMethodT (int i);
  CSEConstants::RenumberMethodT int2RenumberMethodT (int i);
  CSEConstants::SplitMethodT int2SplitMethodT (int i);
};
}
#endif
