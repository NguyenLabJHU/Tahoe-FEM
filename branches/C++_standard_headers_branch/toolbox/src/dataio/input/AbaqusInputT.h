/* $Id: AbaqusInputT.h,v 1.15 2003-11-10 22:14:22 cjkimme Exp $ */
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

namespace Tahoe {

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

	/** ids for all nodes. ids for all nodes in the coordinate list. ids may not be
	 * compact or ordered */
	virtual void ReadNodeID(iArrayT& node_id); 

	/** read coordinates */
	virtual void ReadCoordinates(dArray2DT& coords);
	
	/** read coordinates and ids. ids for all nodes in the coordinate list. 
	 * ids may not be compact or ordered */
	virtual void ReadCoordinates(dArray2DT& coords, iArrayT& node_id);

  int  NumGlobalElements (void) const;
  int  NumElements (const StringT& name);
  int  NumElementNodes (const StringT& name);
  int  NumElementQuadPoints (const StringT& name);
  void ReadAllElementMap (iArrayT& elemmap);
  void ReadGlobalElementMap (const StringT& name, iArrayT& elemmap);
  void ReadGlobalElementSet (const StringT& name, iArrayT& map);
  void ReadConnectivity (const StringT& name, iArray2DT& connects);
  void ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode);

  int  NumNodesInSet (const StringT& name);
  void ReadNodeSet (const StringT& name, iArrayT& nodes);

  bool AreSideSetsLocal (void) const;
  int  NumSidesInSet (const StringT& setname) const;
  StringT SideSetGroupName (const StringT& setname) const;
  void ReadSideSetLocal (const StringT& setname, iArray2DT& sides) const;
  void ReadSideSetGlobal (const StringT& setname, iArray2DT& sides) const;

  void QARecords (ArrayT<StringT>& records);

  int  NumTimeSteps (void) const;
  void ReadTimeSteps (dArrayT& steps);

  int  NumNodeVariables (void) const;
  int  NumElementVariables (void) const;
  int  NumQuadratureVariables (void) const;

  void ReadNodeLabels (ArrayT<StringT>& nlabels) const;
  void ReadElementLabels (ArrayT<StringT>& elabels) const;
  void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;  

  void NodeVariablesUsed (const StringT& name, iArrayT& used);
  void ElementVariablesUsed (const StringT& name, iArrayT& used);
  void QuadratureVariablesUsed (const StringT& name, iArrayT& used);  

  void ReadAllNodeVariable (int step, int varindex, dArrayT& values);
  void ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values);
  void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  void ReadNodeVariables (int step, const StringT& elsetname, dArray2DT& nvalues);
  void ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& nvalues);

  void ReadAllElementVariable (int step, int varindex, dArrayT& values);
  void ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values);
  void ReadAllElementVariables (int step, dArray2DT& evalues);
  void ReadElementVariables (int step, const StringT& name, dArray2DT& evalues);

  void ReadAllQuadratureVariable (int step, int varindex, dArrayT& values);
  void ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values);
  void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  void ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues);

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
inline  int  AbaqusInputT::NumSidesInSet (const StringT& setname)  const
{
#pragma unused (setname)
  return 0; 
}
inline  StringT AbaqusInputT::SideSetGroupName (const StringT& setname)  const
{ 
#pragma unused (setname)
  StringT name ("");
  return name; 
}
inline  void AbaqusInputT::ReadSideSetLocal (const StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free();
}
inline  void AbaqusInputT::ReadSideSetGlobal (const StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free();
}
inline int AbaqusInputT::NumDimensions (void) const { return 3; }
inline int AbaqusInputT::NumElements (const StringT& name) { return fData.NumElements (name); }
inline int AbaqusInputT::NumElementNodes (const StringT& name) { return fData.NumElementNodes (name); }
inline int AbaqusInputT::NumElementQuadPoints (const StringT& name) { return fData.NumElementQuadPoints (name); }
inline int AbaqusInputT::NumQuadratureVariables (void) const { return fData.NumQuadratureVariables (); }
inline int AbaqusInputT::NumNodeVariables (void) const { return fData.NumNodeVariables (); }
inline int AbaqusInputT::NumElementVariables (void) const { return fData.NumElementVariables (); }
inline int AbaqusInputT::NumTimeSteps (void) const { return (fNumModes > 0) ? fNumModes : fNumTimeSteps; }
inline int AbaqusInputT::NumNodesInSet (const StringT& name) { return fData.NumNodesInSet (name); }
inline void AbaqusInputT::ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode) { fData.GeometryCode (name, geocode); }
inline void AbaqusInputT::QARecords (ArrayT<StringT>& records) { fData.VersionNotes (records); }

} // namespace Tahoe 
#endif
