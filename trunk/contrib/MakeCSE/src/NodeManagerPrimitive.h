// file: NodeManagerPrimitive.h

// created: SAW 10/07/99

#ifndef _NODEMANAGERPRIMITIVE_H_
#define _NODEMANAGERPRIMITIVE_H_

#include "dArray2DT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"

class MakeCSEFEManager;
class GlobalEdgeFinderT;
class MakeCSEIOManager;

using namespace Tahoe;

class NodeManagerPrimitive
{
 public:
  NodeManagerPrimitive (ostream& out, int comments, MakeCSEFEManager& FEM);

  void Initialize (MakeCSEIOManager& input);

  /* returns new node number */
  int AddCoord (const dArrayT& p);
  int AddDuplicateCoord (const int oldnode);

  void AddNodeSet (int setID, const ArrayT<int>& nodes, int transfermethod);

  /* returns old node number, given the new node number */
  int OriginalNode (const int node) const;

  void MapNodeSets (const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E);

  void Renumber (int option, iArrayT& map);

  /* accessors */
  int NumNodes (void) const;
  const dArray2DT& InitialCoordinates (void) const;
  int NumNodeSets (void) const;
  iArrayT& NodeSet (int set) const;
  int NodeSetID (int set) const;

  void RegisterOutput (MakeCSEIOManager& theIO);
  
  enum TransferMethods { kSurface1 = 0,
			 kSurface2,
			 kMap,
			 kSplit };

 private:
  void EchoCoordinates (MakeCSEIOManager& theInput);
  void EchoNodeSets (MakeCSEIOManager& theInput);

  void SurfaceNodeSet (iArrayT& set, bool surface1, const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E);
  void MapNodeSet (iArrayT& set, const ArrayT<int>& surface1facets, GlobalEdgeFinderT &E);
  void Split (iArrayT& set);

  void RemoveRepeats (ArrayT<int>& n) const;

 private:
  ostream& out;
  int fPrintUpdate;
  MakeCSEFEManager* theBoss;

  dArray2DT fCoordinates;
  int fNumInitCoordinates;

  /* MakeCSE items */
  ArrayT<iArrayT> fNodeSetData;
  iArrayT fNodeSetID;
  iArrayT fNew2Old;
  iAutoArrayT fSplitNodes;
  ArrayT<iAutoArrayT> fOld2New;

  /* method of node group transfer */
  iArrayT fTransMethods;
};

inline int NodeManagerPrimitive::NumNodes (void) const { return fCoordinates.MajorDim(); }
inline const dArray2DT& NodeManagerPrimitive::InitialCoordinates (void) const { return fCoordinates; }
inline int NodeManagerPrimitive::NumNodeSets (void) const { return fNodeSetData.Length(); }
inline iArrayT& NodeManagerPrimitive::NodeSet (int set) const { return fNodeSetData[set]; }
inline int NodeManagerPrimitive::NodeSetID (int set) const { return fNodeSetID[set]; }

#endif
