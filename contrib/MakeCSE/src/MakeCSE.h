/* $Id: MakeCSE.h,v 1.4 2002-10-08 20:51:50 paklein Exp $ */
/* created      : SAW (04/22/99) */
#ifndef _MAKECSEPRIMITIVE_H_
#define _MAKECSEPRIMITIVE_H_

#include "ClockT.h"
#include "NodeManagerPrimitive.h"
#include "GlobalEdgeFinderT.h"
#include "MakeCSE_ElementBaseT.h"

class MakeCSE_FEManager;

using namespace Tahoe;

/** This class uses the regular mesh to put CSEs into the specified
 * element block, tacks the created nodes on the end of the coordinate list. */
class MakeCSE
{
 public:

  enum ZoneEdgeType { kSingleZE = 1,
		      kDoubleZE = 2,
		      kMixSingZE = 3,
		      kMixDoubZE = 4 };

  MakeCSE (ostream& log, GlobalEdgeFinderT& Edger);
  ~MakeCSE (void);

  void Initialize (MakeCSE_IOManager& theInput, MakeCSE_FEManager& FEM, int comments);
  void Create (void);

 private:

  // gather and initialze data
  void SetFE (MakeCSE_FEManager& FEM);
  void SetInput (MakeCSE_IOManager& theInput);
  void InitializeContact (MakeCSE_IOManager& theInput);
  void CollectFacets (MakeCSE_IOManager& theInput, const iArrayT& facetdata);
  void CollectSingleNodes (MakeCSE_IOManager& theInput);
  void CollectZones (MakeCSE_IOManager& theInput, const iArrayT& zonedata);
  void CollectBoundaries (const iArrayT& boundarydata);

  void InitializeSet (int group, int CSEgroup);
  void AddElements (int num, int CSEgroup);
  void InitializeFacet (int elem, int face, int group, int cse, int csgroup);
  void RemoveSingleNodes (const ArrayT<int>& singles);

  void RenumberFaceNodes (void);
  void MakeList (int node, iAutoArrayT& elems, iAutoArrayT& faces, iAutoArrayT& hit_elems);
  void FindNeighbors (int elm, int face, int node, iAutoArrayT& fElements, iAutoArrayT& fFaces);
  void CheckNeighbor (int elocal, int face, int node, iAutoArrayT& nelems, iAutoArrayT& nfaces);
  void ReduceList (iAutoArrayT& hit_elems, const iAutoArrayT& elems, const iArrayT& checkelems) const;
  int ReNumber (int node, const ArrayT<int>& fElements, const ArrayT<int>& fFaces);

  void CollectMassLessNodes (void);
  void CollectSurfaceData (void);
  void RemoveRepeats (ArrayT<int>& ints) const;

  void PrintControlEnd (ostream& o);

 private: 
  ostream&              out;
  bool                  fPrintUpdate;

  // sides to insert along (elem, facet)
  iAutoArrayT           fSurface1Facets;
  iAutoArrayT           fSurface2Facets;

  // links to data
  GlobalEdgeFinderT*    theEdger;
  NodeManagerPrimitive* theNodes;
  ArrayT<MakeCSE_ElementBaseT*> theElements;
  ArrayT<MakeCSE_ElementBaseT*> theCSEs;
  ClockT                cClock;

  // initial number of nodes and elements before adding CSEs
  int                   fNumStartElements;
  int                   fNumStartNodes;

  // list of nodes to examine for splitting
  iAutoArrayT           fPotentialSplitNodes;

  // courtesy output for user
  iAutoArrayT           fNoMassNodes;

  // list of CSE block ID's to prep for contact surfaces 
  iArrayT               fContact;

  // current id value for adding a set
  int                   fNSetID;
  int                   fSSetID;
};

#endif
