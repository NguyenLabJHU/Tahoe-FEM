// $Id: MakeCSE_ElementBaseT.h,v 1.7 2003-09-10 21:31:08 paklein Exp $
// created: SAW 10/06/00
#ifndef _MakeCSE_ELEMENTBASET_H_
#define _MakeCSE_ELEMENTBASET_H_

#include "iArray2DT.h"
#include "StringT.h"
#include "GeometryT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "sArrayT.h"

namespace Tahoe {

class MakeCSE_IOManager;
class OutputBaseT;
class ModelManagerT;
template <class TYPE> class AutoArrayT;

/** Class similar to ElementBaseT, however, it also has facet data capabilities */
class MakeCSE_ElementBaseT
{
 public:
  /** Constructor */
  MakeCSE_ElementBaseT (ostream& fMainOut, const StringT& ID);

  /** destructor */
  virtual ~MakeCSE_ElementBaseT (void);

  /** read from input that group's connectivity and set facet data */
  void Initialize (ModelManagerT& model, MakeCSE_IOManager& theInput);

  /** no connectivity data exists, element node numbers will be appended later,
      facet data is initialized */
  virtual void Initialize (GeometryT::CodeT geocode, int numnodes);

  /** add elements, reallocates space and initializes new space to kNotSet */
  void AddElements (int numelems);

  /** set the node numbers for an element
   * \param e1local element number, local numbering
   * \param nodes element node numbers, global numbering */
  virtual void SetNodes (int e1local, const iArrayT& nodes);

  /** return list of faces using the node specified
   * \param e1local element number, local numbering
   * \param node node index, global numbering
   * \param faces facet indexes, local numbering */
  void FacesWithNode (int e1local, int node, iArrayT& faces) const;

  /** returns true if element face uses the node specified 
   * \param e1local element index, local numbering
   * \param f1 facet index, local numbering
   * \param node node index, global numbering 
   * \return true if facet has node */
  bool FaceHasNode (int e1local, int f1, int node) const;

  /** replace one node number one a facet of an element 
   * \param e1local element index, local numbering
   * \param f1 facet index, local numbering
   * \param oldnode old node index, global numbering 
   * \param newnode new node index, global numbering */
  void ResetOneFaceNode (int e1global, int f1, int oldnode, int newnode);

  /** add a side set to the list, or append facets to preexisting set
   * \param setID either an existing set ID or new ID
   * \param sides side set facets, globally numbered */
  void AddSideSet (const StringT& setID, const iArray2DT& sides);

  /** renumber element node numbers using map */
  void Renumber (const iArrayT& map);

  // accessors
  int NumElements (void) const; /**< number of elements, increasing as CSE's are created */
  int NumElemFaces (void) const; /**< number of element facets */
  virtual void CSElemFaces (iArrayT& faces) const; /**< returns CSE facets if class is MakeCSE_CSEBaseT, otherwise nothing happens */
  GeometryT::CodeT GeometryCode (void) const; /**< geometry code */
  const StringT& GroupNumber (void) const; /**< group id */
  int NumSideSets (void) const; /**< number of side sets */
  const StringT& SideSetID (int index) const; /**< return side set id */

  /** number of face nodes, accounts for elements with differing number 
     of face nodes */
  int NumFaceNodes (int face) const;

  /** returns element nodes for one element
   * \param e1local element index, local numbering
   * \param nodes element nodes returned */
  void ElementNodes (int e1local, iArrayT& nodes) const; 

  /** returns element nodes for one element
   * \param e1 element index, local numbering 
   * \param f1 facet indes, local numbering
   * \param nodes facet node numbers */
  void FaceNodes (int e1, int f1, iArrayT& nodes) const; 

  /** same as FaceNodes, however, for higher element types, only 
      vertex nodes are returned */
  void AbbrFaceNodes (int e1, int f1, iArrayT& nodes) const;

  /** returns globally numbered nodes used by the element group */
  void NodesUsed (iArrayT& nodes) const;

  /** send output sets and side sets to output manager */
  void RegisterOutput (OutputBaseT& output);

  /** returns true if side set is contained within this element group */
  bool CheckSideSet (const iArray2DT& sides) const;

  /** checks validity range, prints error message, returns true/false */
  bool IsElementValid (int e1local) const;

  /** checks validity range, prints error message, returns true/false */
  bool IsFaceValid (int face) const;

 protected:
  /** store connectivity locally, derived classes are allowed manipulations here */
  virtual void EchoConnectivity (ModelManagerT& theInput);

  /** read connectivity data from model manager */
  void ReadConnectivity (ModelManagerT& theInput, GeometryT::CodeT& geocode, iArray2DT& conn) const;

  /** after storing data, set up facet data */
  void InitializeConnectivity (void);

  /** store side set data locally, derived classes are allowed manipulations here */
  virtual void EchoSideSets (ModelManagerT& model, MakeCSE_IOManager& theInput);

  /** read side set data from model manager */
  void ReadSideSetData (ModelManagerT& model, MakeCSE_IOManager& theInput, ArrayT<iArray2DT>& sides);

  /** verify all side sets are contained within the element group connectivity */
  void CheckAllSideSets (void);

  /** determines facenode map from GeometryT */
  void SetFace (void);

  /** echo data to log */
  void PrintControlData (void) const;

 protected: // share with derived classes
  ostream&          out; /**< log */
  const StringT     fGroupID; /**< group id */
  iArray2DT         fNodeNums; /**< connectivity */
  int               fNumElemNodes; /**< number of nodes per element */
  GeometryT::CodeT  fGeometryCode; /**< element geometry code */
  ArrayT<iArray2DT> fSideSetData; /**< side set data, globally numbered */
  sArrayT           fSideSetID; /**< side set ids */

 private:
  // data from shape class
  ArrayT<iArrayT>    fFacetNodes; /**< for each facet, map of nodes on that facet */
  ArrayT<iAutoArrayT> fRevFacetNodes; /**< reverse of fFacetNodes */
  ArrayT<iArrayT>     fVertexFaceNodes; /**< same as fFacetNodes, except truncated to only vertex nodes for higher order elements */

};

inline int MakeCSE_ElementBaseT::NumElements (void) const { return fNodeNums.MajorDim(); }
inline int MakeCSE_ElementBaseT::NumElemFaces (void) const { return fFacetNodes.Length(); }
inline void MakeCSE_ElementBaseT::CSElemFaces (iArrayT&) const { };
inline GeometryT::CodeT MakeCSE_ElementBaseT::GeometryCode (void) const { return fGeometryCode; }
inline const StringT& MakeCSE_ElementBaseT::GroupNumber (void) const { return fGroupID; }
inline int MakeCSE_ElementBaseT::NumSideSets (void) const { return fSideSetID.Length(); }
inline const StringT& MakeCSE_ElementBaseT::SideSetID (int index) const { return fSideSetID[index]; } 

inline int MakeCSE_ElementBaseT::NumFaceNodes (int face) const 
{ 
  if(IsFaceValid (face)) 
    return fFacetNodes[face].Length(); 
  return -1;
}
}

#endif
