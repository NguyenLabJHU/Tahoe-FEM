#ifndef _PARAMETERFILEMANAGER_T_H_
#define _PARAMETERFILEMANAGER_T_H_

#include "MakeCSE_IOManager.h"

namespace Tahoe {

class ParameterFileManagerT : public MakeCSE_IOManager
{
 public:
  ParameterFileManagerT (const StringT& infile);
  void Initialize (void); /**< verify file existance */

  void InputFormat (IOBaseT::FileTypeT &format, StringT& name); /**< set and echo */
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
  bool AdvanceTo (ifstreamT& in, const StringT& key) const; /**< advance in file to keyword */
  void ReadIDValues (ifstreamT& in, sArrayT& names) const; /**< generic read */
  void ReadID_Parameter (ifstreamT& in, sArrayT& names, iArrayT& params) const; /**< generic read */
  void CheckIDList (sArrayT& names, int numcols, int check) const;
  void CheckIDList (sArrayT& names, iArrayT& itemp) const;

 private:
  const StringT fInFile;
};
}
#endif
