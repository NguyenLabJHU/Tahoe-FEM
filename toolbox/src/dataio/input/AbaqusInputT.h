/* $Id: AbaqusInputT.h,v 1.10 2002-01-23 20:01:58 sawimme Exp $ */
/* created: sawimme (05/18/1998) */

#ifndef _ABAQUSINPUT_T_H_
#define _ABAQUSINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "AbaqusResultsT.h"
#include "StringT.h"
#include "iArray2DT.h"
#include "iArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"

/** If SetLabelName does not have a case for the variable being read, 
 * a default variable name is used. As new variables are added to 
 * AbaqusT::VariableKey, they should also be added to SetLabelName */
class AbaqusInputT : public InputBaseT
{
 public:
  AbaqusInputT (ostream& out);

  bool Open (const StringT& file);
  void Close (void);

  void ElementGroupNames (ArrayT<StringT>& groupnames) const;
  void NodeSetNames (ArrayT<StringT>& nodenames) const;
  void SideSetNames (ArrayT<StringT>& sidenames) const;

  int  NumElementGroups (void) const;
  int  NumSideSets (void) const;
  int  NumNodeSets (void) const;

  int  NumNodes (void) const;
  int  NumDimensions (void) const;
  void ReadNodeMap (iArrayT& nodemap);
  void ReadCoordinates (dArray2DT& coords);
  void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  int  NumGlobalElements (void) const;
  int  NumElements (StringT& name);
  int  NumElementNodes (StringT& name);
  int  NumElementQuadPoints (StringT& name);
  void ReadAllElementMap (iArrayT& elemmap);
  void ReadGlobalElementMap (StringT& name, iArrayT& elemmap);
  void ReadGlobalElementSet (StringT& name, iArrayT& map);
  void ReadConnectivity (StringT& name, iArray2DT& connects);
  void ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode);

  int  NumNodesInSet (StringT& name);
  void ReadNodeSet (StringT& name, iArrayT& nodes);

  bool AreSideSetsLocal (void) const;
  int  NumSidesInSet (StringT& setname) const;
  StringT SideSetGroupName (StringT& setname) const;
  void ReadSideSetLocal (StringT& setname, iArray2DT& sides) const;
  void ReadSideSetGlobal (StringT& setname, iArray2DT& sides) const;

  void QARecords (ArrayT<StringT>& records);

  int  NumTimeSteps (void) const;
  void ReadTimeSteps (dArrayT& steps);

  int  NumNodeVariables (void) const;
  int  NumElementVariables (void) const;
  int  NumQuadratureVariables (void) const;

  void ReadNodeLabels (ArrayT<StringT>& nlabels) const;
  void ReadElementLabels (ArrayT<StringT>& elabels) const;
  void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;  

  void NodeVariablesUsed (StringT& name, iArrayT& used);
  void ElementVariablesUsed (StringT& name, iArrayT& used);
  void QuadratureVariablesUsed (StringT& name, iArrayT& used);  

  void ReadAllNodeVariable (int step, int varindex, dArrayT& values);
  void ReadNodeVariable (int step, StringT& name, int varindex, dArrayT& values);
  void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  void ReadNodeVariables (int step, StringT& elsetname, dArray2DT& nvalues);
  void ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues);

  void ReadAllElementVariable (int step, int varindex, dArrayT& values);
  void ReadElementVariable (int step, StringT& name, int varindex, dArrayT& values);
  void ReadAllElementVariables (int step, dArray2DT& evalues);
  void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);

  void ReadAllQuadratureVariable (int step, int varindex, dArrayT& values);
  void ReadQuadratureVariable (int step, StringT& name, int varindex, dArrayT& values);
  void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues);

private:
  void SetLabelName (const iArrayT& key, const iArrayT& dims, ArrayT<StringT>& name) const;
  void MapOffset (ArrayT<int>& set, const iArrayT& map) const;
  void NodesUsed (const nArrayT<int>& connects, iArrayT& nodesused) const;

 private:
  AbaqusResultsT fData;

  int fNumElements;
  int fNumNodes;
  int fNumTimeSteps;
  int fNumModes;
};

inline void AbaqusInputT::SideSetNames (ArrayT<StringT>& names) const { names.Free (); }
inline int AbaqusInputT::NumElementGroups (void) const { return fData.NumElementSets(); }
inline int AbaqusInputT::NumSideSets (void) const { return 0; }
inline int AbaqusInputT::NumNodeSets (void) const { return fData.NumNodeSets (); }
inline int AbaqusInputT::NumNodes (void) const { return fNumNodes; }
inline int AbaqusInputT::NumGlobalElements (void) const { return fNumElements; }
inline bool AbaqusInputT::AreSideSetsLocal (void) const { return false; }
inline  int  AbaqusInputT::NumSidesInSet (StringT& setname)  const
{ 
#pragma unused (setname)
  return 0; 
}
inline  StringT AbaqusInputT::SideSetGroupName (StringT& setname)  const
{ 
#pragma unused (setname)
  StringT name ("");
  return name; 
}
inline  void AbaqusInputT::ReadSideSetLocal (StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free();
}
inline  void AbaqusInputT::ReadSideSetGlobal (StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free();
}
inline int AbaqusInputT::NumDimensions (void) const { return 3; }
inline int AbaqusInputT::NumElements (StringT& name) { return fData.NumElements (name); }
inline int AbaqusInputT::NumElementNodes (StringT& name) { return fData.NumElementNodes (name); }
inline int AbaqusInputT::NumElementQuadPoints (StringT& name) { return fData.NumElementQuadPoints (name); }
inline int AbaqusInputT::NumQuadratureVariables (void) const { return fData.NumQuadratureVariables (); }
inline int AbaqusInputT::NumNodeVariables (void) const { return fData.NumNodeVariables (); }
inline int AbaqusInputT::NumElementVariables (void) const { return fData.NumElementVariables (); }
inline int AbaqusInputT::NumTimeSteps (void) const { return (fNumModes > 0) ? fNumModes : fNumTimeSteps; }
inline int AbaqusInputT::NumNodesInSet (StringT& name) { return fData.NumNodesInSet (name); }
inline void AbaqusInputT::ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode) { fData.GeometryCode (name, geocode); }
inline void AbaqusInputT::QARecords (ArrayT<StringT>& records) { fData.VersionNotes (records); }

#endif
