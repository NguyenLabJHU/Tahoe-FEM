/* $Id: ModelManagerT.h,v 1.4.2.4 2001-10-11 19:58:19 sawimme Exp $ */
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

  /* Read format and file name from file, casts and initializes InputBaseT */
  void Initialize (ifstreamT& in);

  /* give file name and format, casts and initializes InputBaseT */
  void Initialize (const IOBaseT::FileTypeT format, const StringT& database);

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
     allows code to read from input file but not care about input format */
  void SideSetList (ifstreamT& in, iArrayT& indexes);

  /* read IC/KBC/FBC card type data
     allows code to read from input file but not care about input format */
  int ModelManagerT::ReadCards (ifstreamT& in, ostream& out, ArrayT<iArrayT>& nodes, iArray2DT& data, dArrayT& value);

  /* access */
  void CoordinateDimensions (int& length, int& dof) const;
  const dArray2DT& CoordinateReference (void) const; /* return coords */
  const dArray2DT& Coordinates (void); /* reads coordinates if not yet read and returns coords */
  void ReadCoordinates (void); /* reads coordinates, assumes coords reference is set separately */

  int NumElementGroups (void) const;
  void ElementGroupNames (ArrayT<StringT>& names) const;
  int ElementGroupIndex (const StringT& name) const;
  void ElementGroupDimensions (int index, int& numelems, int& numelemnodes) const;
  GeometryT::CodeT ElementGroupGeometry (int index) const;
  const iArray2DT& ElementGroup (int index);

  void AllNodeMap (iArrayT& map);
  void AllElementMap (iArrayT& map);
  void ElementMap (StringT& name, iArrayT& map);

  int NumNodeSets (void) const;
  void NodeSetNames (ArrayT<StringT>& names) const;
  int NodeSetIndex (const StringT& name) const;
  int NodeSetLength (int index) const;
  const iArrayT& NodeSet (int index);
  void ManyNodeSets (const iArrayT& indexes, iArrayT& nodes);

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

  /* This closes link to InputBaseT, it does not clear data */
  void CloseModel (void);

  /* checks to see if external file or actual data */
  ifstreamT& OpenExternal (ifstreamT& in, ifstreamT& in2, ostream& out, bool verbose, const char* fail) const;
  
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

inline int ModelManagerT::NumElementGroups (void) const { return fNumElementSets; }
inline int ModelManagerT::NumNodeSets (void) const { return fNumNodeSets; }
inline int ModelManagerT::NumSideSets (void) const { return fNumSideSets; }
#endif
