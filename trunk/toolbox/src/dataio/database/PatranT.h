/* $Id: PatranT.h,v 1.10 2002-07-23 11:39:12 sawimme Exp $ */
/* created: sawimme (05/17/2001)  */

#ifndef _PATRAN_T_H_
#define _PATRAN_T_H_

/* direct members */
#include "StringT.h"
#include "iArray2DT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class iArrayT;
class dArray2DT;
class dArrayT;

class PatranT
{
 public:
  enum NamedTypes { kNCPoint = 1, kNCCurve = 2,
		    kNCPatch = 3, kNCHyperPatch = 4,
		    kNCNode = 5, 
		    kNCLine = 6, kNCLine2 = 106, kNCLine3 = 206,
		    kNCTriangle = 7, kNCTriangle2 = 107, kNCTriangle3 = 207,
		    kNCQuad = 8, kNCQuad2 = 108, kNCQuad3 = 208,
		    kNCTet = 9, kNCTet2 = 109, kNCTet3 = 209,
		    kNCWedge = 11, kNCWedge2 = 111, kNCWedge3 = 211,
		    kNCHex = 12, kNCHex2 = 112, kNCHex3 = 212 };

  enum ElementTypes { kLine = 2,
		      kTriangle = 3,
		      kQuadrilateral = 4,
		      kTetrahedron = 5,
		      kPentahedron = 7,
		      kHexahedron = 8 };

  PatranT (ostream &messge_out);
  ~PatranT (void);

  bool OpenRead (const StringT& filename);

  /* accessors */
  const StringT& Filename (void) const;
  void VersionNotes (ArrayT<StringT>& records) const;
  int NumNodes (void) const;
  int NumElements (void) const;
  int NumDimensions (void) const;
  int NumNamedComponents (void) const;
  bool NamedComponents (ArrayT<StringT>& names) const;
  bool NumNodesInSet (const StringT& title, int& num) const;
  bool ReadGlobalNodeMap (iArrayT& map) const;
  bool ReadGlobalElementMap (iArrayT& map) const;
  bool ReadCoordinates (dArray2DT& coords, int dof) const;
  bool ReadElementBlockDims (const StringT& title, int& num_elems, int& num_elem_nodes) const;
  bool ReadConnectivity (const StringT& title, int& namedtype, iArray2DT& connects) const;
  bool ReadAllElements (ArrayT<iArrayT>& connects, iArrayT& elementtypes) const;
  bool ReadElementSet (const StringT& title, int& namedtype, iArrayT& elems) const;
  bool ReadElementSetMixed (const StringT& title, iArrayT& namedtype, iArrayT& elems) const;
  bool ReadDistLoadSetDims (int setID, int& num_elems) const;
  bool ReadDistLoadSet (int setID, iArray2DT& facets) const;
  bool ReadNodeSet (const StringT& title, iArrayT& nodes) const;
  bool ReadNodeSets (const ArrayT<StringT>& title, iArrayT& nodes) const;

  /* write geometry file */
  bool WriteHeader (ostream& out, int numnodes, int numelems, const StringT& title) const;
  bool WriteCoordinates (ostream& out, const dArray2DT& coords, int firstnodeID) const;
  bool WriteElements (ostream& out, const iArray2DT& elems, const iArrayT& elemtypes, int firstelemID) const;
  bool WriteNamedComponent (ostream& out, const StringT& name, int ID, const iArray2DT& comps) const;
  bool WriteGeometryPoints (ostream& out, const dArray2DT& points, int firstptiD) const;
  bool WritePairPointCurve (ostream& out, int curveID, int ID1, int ID2, const dArrayT& coord1, const dArrayT& coord2) const;
  bool WriteClosure (ostream& out) const;

  /* write results data files */
  bool WriteNodalVariables (const StringT& filename, const iArrayT& ids, const dArray2DT& values, const ArrayT<StringT>& titles) const;
  bool WriteElementVariables (const StringT& filename, const iArrayT& ids, const ArrayT<PatranT::ElementTypes>& shapes, const dArray2DT& values, const ArrayT<StringT>& titles) const; 

 private:
  enum PacketT { kTitle = 25,
		 kSummary = 26,
		 kNode = 1,
		 kElement = 2,
		 kMaterial = 3,
		 kElementProps = 4,
		 kFrame = 5,
		 kDistLoads = 6,
		 kNodeForce = 7,
		 kNodeDisp = 8,
		 kNamedComponents = 21,
		 kGridData = 31,
		 kLineData = 32 };

  enum ShapeT { kBarShape = 2,
		kTriShape = 3,
		kQuadShape = 4,
		kTetShape = 5,
		kWedgeShape = 7,
		kHexShape = 8 };

  enum NTypeT { kNodeType = 5,
		kBarType = 6,
		kTriType = 7,
		kQuadType = 8,
		kTetType = 9,
		kWedgeType = 11, 
		kHexType = 12,
		kBarTypeA = 106,
		kTriTypeA = 107,
		kQuadTypeA = 108,
		kTetTypeA = 109,
		kWedgeTypeA = 111, 
		kHexTypeA = 112,
		kBarTypeB = 206,
		kTriTypeB = 207,
		kQuadTypeB = 208,
		kTetTypeB = 209,
		kWedgeTypeB = 211, 
		kHexTypeB = 212 };

  enum GeometryFormatT { hwidth = 8, cwidth = 16, prec = 9};
  enum ResultsFormatT { kString = 80, kHeaderInt = 9, kDataInt = 8,
			kHeaderDouble = 15, kDataDouble=13,
			kHeaderPrec = 6, kDataPrec = 7 };

  void ScanFile (void);
  int LocateNamedComponent (const StringT &title) const;
  bool AdvanceTo (ifstream& in, int target, int &ID, int &IV, int &KC) const;
  void ClearPackets (ifstream &in, int KC) const;

  bool WritePacketHeader (ostream& out, int tag, int ID, int IV, int KC, iArrayT n) const;

 private:
  StringT file_name;
  ostream &fMessage;

  ArrayT<StringT> fNamedComponents; /**< Named Components */
  ArrayT<iArray2DT> fNamedComponentsData; /**< Named Component Data: Type, ID */
};

inline const StringT& PatranT::Filename (void) const { return file_name; }

} // namespace Tahoe 
#endif
