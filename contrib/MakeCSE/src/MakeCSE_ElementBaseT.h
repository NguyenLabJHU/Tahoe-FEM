// file: ElementBaseT.h

// created: SAW 10/06/00

#ifndef _MakeCSE_ELEMENTBASET_H_
#define _MakeCSE_ELEMENTBASET_H_

#include "iArray2DT.h"
#include "StringT.h"
#include "GeometryT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "AutoArrayT.h"
#include "ModelManagerT.h"
#include "sArrayT.h"

namespace Tahoe {

class MakeCSE_IOManager;

class MakeCSE_ElementBaseT
{
 public:
  MakeCSE_ElementBaseT (ostream& fMainOut, const StringT& ID);

  virtual ~MakeCSE_ElementBaseT (void);

  // read from input that group's connectivity, create element cards
  void Initialize (ModelManagerT& model, MakeCSE_IOManager& theInput);

  // initialize element group for newly created data
  virtual void Initialize (GeometryT::CodeT geocode, int numnodes);

  // add elements, reallocates space and initializes new space to kNotSet
  void AddElements (int numelems);

  // set the node numbers for an element
  virtual void SetNodes (int e1local, const iArrayT& nodes);

  // return list of faces using the node specified
  void FacesWithNode (int e1local, int node, iArrayT& faces) const;

  // returns true if element face uses the node specified
  bool FaceHasNode (int e1local, int f1, int node) const;

  // replace one node number one a facet of an element
  void ResetOneFaceNode (int e1global, int f1, int oldnode, int newnode);

  // add a side set to the list, or append facets to preexisting set
  void AddSideSet (const StringT& setID, const iArray2DT& sides);

  // renumber node numbering
  void Renumber (const iArrayT& map);

  // accessors
  int NumElements (void) const;
  int NumElemFaces (void) const;
  virtual void CSElemFaces (iArrayT& faces) const;
  GeometryT::CodeT GeometryCode (void) const;
  const StringT& GroupNumber (void) const; // exterior numbering
  int NumSideSets (void) const;

  int NumFaceNodes (int face) const;
  void ElementNodes (int e1local, iArrayT& nodes) const;
  void FaceNodes (int e1, int f1, iArrayT& nodes) const;
  void AbbrFaceNodes (int e1, int f1, iArrayT& nodes) const;
  iArray2DT& SideSet (int g) const;
  const StringT& SideSetID (int g) const;
  iArray2DT& Neighbors (void);

  /* accessor for output manager */
  const iArray2DT& Connectivity (void) const;
  void Connectivities (AutoArrayT<const iArray2DT*>& conn, iAutoArrayT& geocodes, iAutoArrayT& change);

  void NodesUsed (iArrayT& nodes) const;
  void RegisterOutput (MakeCSE_IOManager& theIO);
  void WriteOutput (void) const;

  /* returns true if side set is contained within this element group */
  bool CheckSideSet (const iArray2DT& sides) const;

  // checks validity range, prints error message, returns true/false
  bool IsElementValid (int e1local) const;
  bool IsFaceValid (int face) const;

 protected:
  virtual void EchoConnectivity (ModelManagerT& theInput);
  void ReadConnectivity (ModelManagerT& theInput, GeometryT::CodeT& geocode, iArray2DT& conn) const;
  void InitializeConnectivity (void);
  virtual void EchoSideSets (ModelManagerT& model, MakeCSE_IOManager& theInput);
  void ReadSideSetData (ModelManagerT& model, MakeCSE_IOManager& theInput, ArrayT<iArray2DT>& sides);
  void CheckAllSideSets (void);

  /* determines facenode map from GeometryT */
  void SetFace (void);

  void PrintControlData (void) const;

 protected: // share with derived classes
  ostream&          out;
  const StringT     fGroupID;    // external numbering
  int               fOutputID;   // output set ID

  iArray2DT         fNodeNums;
  int               fNumElemNodes;
  GeometryT::CodeT  fGeometryCode;
  ArrayT<iArray2DT> fSideSetData;
  sArrayT           fSideSetID;

 private:
  // data from shape class
  ArrayT<iArrayT>    fFacetNodes;
  ArrayT<iAutoArrayT> fRevFacetNodes;
  ArrayT<iArrayT>     fVertexFaceNodes; // only the vertex nodes

};

inline int MakeCSE_ElementBaseT::NumElements (void) const { return fNodeNums.MajorDim(); }
inline int MakeCSE_ElementBaseT::NumElemFaces (void) const { return fFacetNodes.Length(); }
inline void MakeCSE_ElementBaseT::CSElemFaces (iArrayT& faces) const { };
inline GeometryT::CodeT MakeCSE_ElementBaseT::GeometryCode (void) const { return fGeometryCode; }
inline const StringT& MakeCSE_ElementBaseT::GroupNumber (void) const { return fGroupID; }
inline int MakeCSE_ElementBaseT::NumSideSets (void) const { return fSideSetData.Length(); }
inline iArray2DT& MakeCSE_ElementBaseT::SideSet (int g) const { return fSideSetData[g]; }
inline const StringT& MakeCSE_ElementBaseT::SideSetID (int g) const { return fSideSetID[g]; }

inline const iArray2DT& MakeCSE_ElementBaseT::Connectivity (void) const {return fNodeNums; }
inline int MakeCSE_ElementBaseT::NumFaceNodes (int face) const 
{ 
  if(IsFaceValid (face)) 
    return fFacetNodes[face].Length(); 
  return -1;
}
}

#endif
