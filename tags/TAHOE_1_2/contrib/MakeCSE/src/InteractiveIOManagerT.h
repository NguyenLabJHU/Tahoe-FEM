#ifndef _INTERACTIVEIOMANAGER_T_H_
#define _INTERACITVEIOMANAGER_T_H_

#include "MakeCSE_IOManager.h"

namespace Tahoe {

class InteractiveIOManagerT : public MakeCSE_IOManager
{
 public:
  InteractiveIOManagerT (void);
  void Initialize (void); /**< Query user about parameters */
  
  void InputFormat (IOBaseT::FileTypeT &format, StringT& name);
  void OutputFormat (IOBaseT::FileTypeT &format, StringT& name);
  bool Verbose (void);
  void Facets (sArrayT& names);
  void Zones (sArrayT& names);
  void Boundaries (sArrayT& names);
  CSEConstants::ZoneEdgeT ZoneMethod (void);
  void ZoneEdgeNodeSets (sArrayT& names);
  void Contact (sArrayT& names);
  void SingleNodes (sArrayT& names);
  void BlockToNode (sArrayT& names);
  void NodeSetsMapped (sArrayT& names, ArrayT<CSEConstants::NodeMapMethodT>& meths);
  void SideSetsMapped (sArrayT& names);
  CSEConstants::RenumberMethodT RenumberMethod (void);
  void SplitBlocks (sArrayT& names, ArrayT<CSEConstants::SplitMethodT>& meths);

 private:
  void Method (void);
  void ReadIDValues (const sArrayT& q, sArrayT& names); /** generic reading of StringT id values */
  void ReadID_Parameter (const sArrayT& q, sArrayT& names, iArrayT& vals); /** read mapped node sets and their methods */

 private:
  ofstreamT fEchoInput; /**< input file to write based on interactive answers */
  bool fEcho; /**< are we creating fEchoInput */

  CSEConstants::CSEMethodT fMethod;
  sArrayT fMethodSets;
};
} //namespace Tahoe
#endif
