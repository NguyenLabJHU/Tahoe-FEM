// file: NodeManagerPrimitive.h

// created: SAW 10/07/99

#ifndef _NODEMANAGERPRIMITIVE_H_
#define _NODEMANAGERPRIMITIVE_H_

#include "dArray2DT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "CSEConstants.h"
#include "ModelManagerT.h"
#include "sArrayT.h"

namespace Tahoe {

class MakeCSE_FEManager;
class GlobalEdgeFinderT;
class MakeCSE_IOManager;
class OutputBaseT;

class NodeManagerPrimitive
{
 public:
  /** constructor */
  NodeManagerPrimitive (ostream& out, bool comments, MakeCSE_FEManager& FEM);

  /** read data from model */
  void Initialize (ModelManagerT& model, MakeCSE_IOManager& input);

  /** add coordinate and return new node number */
  int AddCoord (const dArrayT& p);

  /** add duplicate coordiante and return new node number */
  int AddDuplicateCoord (const int oldnode);

  /** add nodes to set if set exists, otherwise, add set */
  void AddNodeSet (const StringT& setID, const ArrayT<int>& nodes, CSEConstants::NodeMapMethodT transfermethod);

  /** returns old node number, given the new node number (after node has been split or double noded */
  int OriginalNode (const int node) const;

  /** map node sets to new node numbering, after double noding or spliting the nodes */
  void MapNodeSets (const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E);

  /** renumber nodes, so new nodes are inserted into the list */
  void Renumber (CSEConstants::RenumberMethodT option, iArrayT& map);

  /* accessors */
  int NumNodes (void) const; /** current number of nodes */
  const dArray2DT& InitialCoordinates (void) const; /** current coordinate list */
  int NumNodeSets (void) const;
  const StringT& NodeSetID (int index) const;

  /** register coordinates and node sets to output manager
      create node sets of element blocks if desired */
  void RegisterOutput (OutputBaseT& output, MakeCSE_IOManager& theInput);
  
 private:
  /** read coordinates from model */
  void EchoCoordinates (ModelManagerT& theInput);

  /** read node sets from model */
  void EchoNodeSets (ModelManagerT& model, MakeCSE_IOManager& theInput);

  /** map node set using surface method, keep node on either surface 1 or surface 2, where CSE's were inserted */
  void SurfaceNodeSet (iArrayT& set, bool surface1, const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E);

  /** map node set, determine which surface the node is on */
  void MapNodeSet (iArrayT& set, const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E);

  /** map node set, add nodes from both surfaces to the set */
  void Split (iArrayT& set);

  /** remove repeats from node set, could move function to nArrayT? */
  void RemoveRepeats (ArrayT<int>& n) const;

 private:
  ostream& out;
  bool fPrintUpdate;
  MakeCSE_FEManager* theBoss;

  dArray2DT fCoordinates;
  int fNumInitCoordinates;

  /** Temporary output node tags, not currently saved values */
  iArrayT fOutputNodeMap;

  /* MakeCSE items */
  ArrayT<iArrayT> fNodeSetData; /**< node set data */
  sArrayT fNodeSetID; /**< node set ids */

  /* contain index values, not node tags */
  iArrayT fNew2Old; /**< map stating what the original node was before the node was split */
  iAutoArrayT fSplitNodes; /**< nodes that were split */
  ArrayT<iAutoArrayT> fOld2New; /**< reverse of fNew2Old */

  /** method of node group transfer */
  ArrayT<CSEConstants::NodeMapMethodT> fTransMethods;
};

inline int NodeManagerPrimitive::NumNodes (void) const { return fCoordinates.MajorDim(); }
inline const dArray2DT& NodeManagerPrimitive::InitialCoordinates (void) const { return fCoordinates; }
inline  int NodeManagerPrimitive::NumNodeSets (void) const { return fNodeSetID.Length(); }
inline  const StringT& NodeManagerPrimitive::NodeSetID (int index) const { return fNodeSetID[index]; }
}
#endif
