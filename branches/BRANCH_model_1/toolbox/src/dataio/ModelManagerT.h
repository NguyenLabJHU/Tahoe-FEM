/* $Id: ModelManagerT.h,v 1.4.2.7 2001-10-18 21:48:33 sawimme Exp $ */
/* created: sawimme July 2001 */

#ifndef _MODELMANAGER_T_H_
#define _MODELMANAGER_T_H_

/* direct */
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "StringT.h"
#include "GeometryT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
#include "InputBaseT.h"

/* forward */
#include "ios_fwd_decl.h"
class ifstreamT;

class ModelManagerT
{
 public:
  ModelManagerT (ostream& message);
  ~ModelManagerT (void);

  /* Read format and file name from file, casts and initializes InputBaseT if not readonly */
  void Initialize (ifstreamT& in, bool readonly);

  /* give file name and format, casts and initializes InputBaseT */
  void Initialize (const IOBaseT::FileTypeT format, const StringT& database);

  /* echo format and model file to message file */
  void EchoData (ostream& o) const;
  void Format (IOBaseT::FileTypeT& format, StringT& name) const;

  /* Query the user interactively for format and file name, 
     for translator programs,
     casts and initializes InputBaseT */
  void Initialize (void);

  /* registration, arrays to be read later as needed  */
  bool RegisterNodes (int length, int dof);
  bool RegisterElementGroup (const StringT& name, int numelems, int numelemnodes, GeometryT::CodeT code);
  bool RegisterNodeSet (const StringT& name, int length);
  bool RegisterSideSet (const StringT& name, int length, bool local, int groupindex);

  /* register the actual array */
  bool RegisterNodes (dArray2DT& coords);
  bool RegisterElementGroup (const StringT& name, iArray2DT& conn, GeometryT::CodeT code);
  bool RegisterNodeSet (const StringT& name, iArrayT& set);
  bool RegisterSideSet (const StringT& name, iArray2DT& set, bool local, int groupindex);

  /* register the actual array from a file */
  bool RegisterNodes (ifstreamT& in);
  bool RegisterElementGroup (ifstreamT& in, const StringT& name, GeometryT::CodeT code);
  bool RegisterNodeSet (ifstreamT& in, const StringT& name);
  bool RegisterSideSet (ifstreamT& in, const StringT& name, bool local, int groupindex);

  /* reads from input file the coordinate dimensions and coords if kTahoe
     allows code to read from input file but not care about input format */
  void ReadInlineCoordinates (ifstreamT& in);

  /* reads from input file the number of blocks, list of names, and matnums
     if kTahoe, it will read connectivity and register it
     return array of indexes and array of matnums
     allows code to read from input file but not care about input format */
  void ElementBlockList (ifstreamT& in, iArrayT& indexes, iArrayT& matnums);

  /* reads from input file the number of sets and list of names
     if kTahoe, it will read number of nodes and register them
     returns array of indexes
     allows code to read from input file but not care about input format */
  void NodeSetList (ifstreamT& in, iArrayT& indexes);

  /* reads from input file the number of sets and list of names
     if kTahoe, it will read number of facets and register them
     returns array of indexes
     allows code to read from input file but not care about input format 
     multidatabasesets accounts for the inconsistency between node set lists and side sets lists
     many side set lists are only allowed one side set and therefore the num_sets is not read from parameter file */
  void SideSetList (ifstreamT& in, iArrayT& indexes, bool multidatabasesets);

  /* read IC/KBC/FBC card type data
     allows code to read from input file but not care about input format */
  int ReadCards (ifstreamT& in, ostream& out, ArrayT<iArrayT>& nodes, iArray2DT& data, dArrayT& value);

  /* read traction card type data */
  void ModelManagerT::ReadNumTractionLines (ifstreamT& in, int& numlines, int& numsets);
  void ModelManagerT::ReadTractionSetData (ifstreamT& in, int& blockindex, int& setsize);
  void ModelManagerT::ReadTractionSideSet (ifstreamT& in, int& blockindex, iArray2DT& localsides);

  /* access node data */
  void CoordinateDimensions (int& length, int& dof) const;
  const dArray2DT& CoordinateReference (void) const; /* return coords */
  const dArray2DT& Coordinates (void); /* reads coordinates if not yet read and returns coords */
  void ReadCoordinates (void); /* reads coordinates, assumes coords reference is set separately */

  /* use to determine if 3D coordinates with 2D elements
     Patran, Abaqus, EnSight, etc. always store coordinates in 3D */
  bool AreElements2D (void) const;

  /* access element data */
  int NumElementGroups (void) const;
  void ElementGroupNames (ArrayT<StringT>& names) const;
  int ElementGroupIndex (const StringT& name) const;
  void ElementGroupDimensions (int index, int& numelems, int& numelemnodes) const;
  GeometryT::CodeT ElementGroupGeometry (int index) const;
  const iArray2DT& ElementGroup (int index); /* read elements if not yet read and returns elements */

  /* access map data, data is not offset */
  void AllNodeMap (iArrayT& map);
  void AllElementMap (iArrayT& map);
  void ElementMap (StringT& name, iArrayT& map);

  /* access node set data */
  int NumNodeSets (void) const;
  void NodeSetNames (ArrayT<StringT>& names) const;
  int NodeSetIndex (const StringT& name) const;
  int NodeSetLength (int index) const;
  const iArrayT& NodeSet (int index);
  void ManyNodeSets (const iArrayT& indexes, iArrayT& nodes);

  /* access side set data */
  int NumSideSets (void) const;
  void SideSetNames (ArrayT<StringT>& names) const;
  int SideSetIndex (const StringT& name) const;
  int SideSetLength (int index) const;
  const iArray2DT& SideSet (int index) const;
  bool IsSideSetLocal (int index) const;
  int SideSetGroupIndex (int sidesetindex) const;
  void SideSetLocalToGlobal (const int localelemindex, const iArray2DT& local, iArray2DT& global);
  void SideSetGlobalToLocal (int& localelemindex, iArray2DT& local, const iArray2DT& global);

  /* modifiers */
  void AddNodes (const dArray2DT& newcoords, iArrayT& new_node_tags, int& newtotalnumnodes);
  void DuplicateNodes (const iArrayT& nodes, iArrayT& new_node_tags, int& newtotalnumnodes);
  void AdjustCoordinatesto2D (void);

  /* This closes link to InputBaseT, it does not clear data */
  void CloseModel (void);

  /* checks to see if external file or actual data */
  ifstreamT& OpenExternal (ifstreamT& in, ifstreamT& in2, ostream& out, bool verbose, const char* fail) const;

  /* access time steps */
  int NumTimeSteps (void);
  void TimeSteps (dArrayT& steps);
  
  /* access node variables */
  int NumNodeVariables (void);
  void NodeLabels (ArrayT<StringT>& labels);
  void AllNodeVariables (int stepindex, dArray2DT& values); /* for all node points */
  void NodeVariables (int stepindex, StringT& elsetname, dArray2DT& values); /* for nodes in element set */
  void NodeSetVariables (int stepindex, StringT& nsetname, dArray2DT& values); /* for nodes in node set */

  /* access element variables */
  int NumElementVariables (void);
  void ElementLabels (ArrayT<StringT>& labels);
  void AllElementVariables (int stepindex, dArray2DT& values); /* for all elements */
  void ElementVariables (int stepindex, StringT& elsetname, dArray2DT& values); /* for element set */

  /* access quadrature variables */
  int NumElementQuadPoints (StringT& elsetname);
  int NumQuadratureVariables (void);
  void QuadratureLabels (ArrayT<StringT>& labels);
  void AllQuadratureVariables (int stepindex, dArray2DT& values); /* for all elements */
  void QuadratureVariables (int stepindex, StringT& elsetname, dArray2DT& values); /* for element set */

 private:
  void ScanModel (const StringT& database);
  bool ScanElements (void);
  bool ScanNodeSets (void);
  bool ScanSideSets (void);
  bool CheckName (const ArrayT<StringT>& list, const StringT& name, const char *type) const;

 protected:
  ostream& fMessage;
  IOBaseT::FileTypeT fFormat;
  InputBaseT *fInput;
  StringT fInputName;

 private:
  /* dimensional information */
  iArrayT fCoordinateDimensions;
  iAutoArrayT fElementLengths;
  iAutoArrayT fElementNodes;
  iAutoArrayT fNodeSetDimensions;
  iAutoArrayT fSideSetDimensions;
  AutoArrayT<bool> fSideSetIsLocal;
  iAutoArrayT fSideSetGroupIndex;

  /* set parameters */
  AutoArrayT<StringT> fElementNames;
  AutoArrayT<StringT> fNodeSetNames;
  AutoArrayT<StringT> fSideSetNames;
  AutoArrayT<GeometryT::CodeT> fElementCodes;

  /* data */
  int fNumElementSets;
  int fNumNodeSets;
  int fNumSideSets;
  dArray2DT fCoordinates;
  AutoArrayT<iArray2DT> fElementSets;
  AutoArrayT<iArrayT> fNodeSets;
  AutoArrayT<iArray2DT> fSideSets;
};

inline void ModelManagerT::Format (IOBaseT::FileTypeT& format, StringT& name) const
{ 
  format = fFormat;
  name = fInputName;
}

inline int ModelManagerT::NumElementGroups (void) const { return fNumElementSets; }
inline int ModelManagerT::NumNodeSets (void) const { return fNumNodeSets; }
inline int ModelManagerT::NumSideSets (void) const { return fNumSideSets; }

inline int ModelManagerT::NumTimeSteps (void) { return fInput->NumTimeSteps(); }
inline void ModelManagerT::TimeSteps (dArrayT& steps) { fInput->ReadTimeSteps (steps); }

inline int ModelManagerT::NumNodeVariables (void) { return fInput->NumNodeVariables (); }
inline void ModelManagerT::NodeLabels (ArrayT<StringT>& labels) { fInput->ReadNodeLabels (labels); }
inline void ModelManagerT::AllNodeVariables (int stepindex, dArray2DT& values) { fInput->ReadAllNodeVariables (stepindex, values); }
inline void ModelManagerT::NodeVariables (int stepindex, StringT& elsetname, dArray2DT& values) { fInput->ReadNodeVariables (stepindex, elsetname, values); }
inline void ModelManagerT::NodeSetVariables (int stepindex, StringT& nsetname, dArray2DT& values) { fInput->ReadNodeSetVariables (stepindex, nsetname, values); }

inline int ModelManagerT::NumElementVariables (void) { return fInput->NumElementVariables (); }
inline void ModelManagerT::ElementLabels (ArrayT<StringT>& labels) { fInput->ReadElementLabels (labels); }
inline void ModelManagerT::AllElementVariables (int stepindex, dArray2DT& values) { fInput->ReadAllElementVariables (stepindex, values); }
inline void ModelManagerT::ElementVariables (int stepindex, StringT& elsetname, dArray2DT& values) { fInput->ReadElementVariables (stepindex, elsetname, values); }

inline int ModelManagerT::NumElementQuadPoints (StringT& name) { return fInput->NumElementQuadPoints(name); }
inline int ModelManagerT::NumQuadratureVariables (void) { return fInput->NumQuadratureVariables (); }
inline void ModelManagerT::QuadratureLabels (ArrayT<StringT>& labels) { fInput->ReadQuadratureLabels (labels); }
inline void ModelManagerT::AllQuadratureVariables (int stepindex, dArray2DT& values) { fInput->ReadAllQuadratureVariables (stepindex, values); }
inline void ModelManagerT::QuadratureVariables (int stepindex, StringT& elsetname, dArray2DT& values) { fInput->ReadQuadratureVariables (stepindex, elsetname, values); }

#endif
