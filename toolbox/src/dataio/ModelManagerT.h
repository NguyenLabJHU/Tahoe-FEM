/* $Id: ModelManagerT.h,v 1.2 2001-08-07 23:11:52 paklein Exp $ */
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

/** still need to add global node and element maps
 * assume for now that side sets are globally numbered */
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

  /* access */
  void CoordinateDimensions (int& length, int& dof) const;
  const dArray2DT& Coordinates (void);

  int NumElementGroups (void) const;
  void ElementGroupNames (ArrayT<StringT>& names) const;
  int ElementGroupIndex (const StringT& name) const;
  void ElementGroupDimensions (int index, int& numelems, int& numelemnodes) const;
  GeometryT::CodeT ElementGroupGeometry (int index) const;
  const iArray2DT& ElementGroup (int index);

  int NumNodeSets (void) const;
  int NodeSetIndex (const StringT& name) const;
  int NodeSetLength (int index) const;
  const iArrayT& NodeSet (int index);

  int NumSideSets (void) const;
  int SideSetIndex (const StringT& name) const;
  int SideSetLength (int index) const;
  const iArray2DT& SideSet (int index) const;
  bool IsSideSetLocal (int index) const;
  int SideSetGroupIndex (int sidesetindex) const;
  void SideSetLocalToGlobal (const int localelemindex, const iArray2DT& local, iArray2DT& global);
  void SideSetGlobalToLocal (int& localelemindex, iArray2DT& local, const iArray2DT& global);

  void CloseModel (void);
  
 private:
  void ScanModel (const IOBaseT::FileTypeT format, const StringT& database);
  bool ScanElements (void);
  bool ScanNodeSets (void);
  bool ScanSideSets (void);

 private:
  ostream& fMessage;
  InputBaseT *fInput;

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
