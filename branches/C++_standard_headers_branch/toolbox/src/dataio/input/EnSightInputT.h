/* $Id: EnSightInputT.h,v 1.12 2002-07-05 22:26:26 paklein Exp $ */
/* created: sawimme (05/18/1998) */

#ifndef _ENSIGHTINPUT_T_H_
#define _ENSIGHTINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "EnSightT.h"
#include "StringT.h"
#include "iArray2DT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
template <class TYPE> class ArrayT;

class iArray2DT;

class EnSightInputT : public InputBaseT
{
public:
  EnSightInputT (ostream& out, bool binary);

  virtual bool Open (const StringT& file);
  virtual void Close (void);

  /* virtual with InputManager base class */
  virtual void ElementGroupNames (ArrayT<StringT>& groupnames) const;
  virtual void SideSetNames (ArrayT<StringT>& sidenames) const;
  virtual void NodeSetNames (ArrayT<StringT>& nodenames) const;

  virtual int  NumElementGroups (void) const;
  virtual int  NumSideSets (void) const;
  virtual int  NumNodeSets (void) const;

  virtual int  NumNodes (void) const;
  virtual int  NumDimensions (void) const;

	/** ids for all nodes. ids for all nodes in the coordinate list. ids may not be
	 * compact or ordered */
	virtual void ReadNodeID(iArrayT& node_id); 

	/** read coordinates */
	virtual void ReadCoordinates(dArray2DT& coords);
	
	/** read coordinates and ids. ids for all nodes in the coordinate list. 
	 * ids may not be compact or ordered */
	virtual void ReadCoordinates(dArray2DT& coords, iArrayT& node_id);

  virtual int  NumGlobalElements (void) const;
  virtual int  NumElements (const StringT& name);
  virtual int  NumElementNodes (const StringT& name);
  virtual int  NumElementQuadPoints (const StringT& name);
  virtual void ReadAllElementMap (iArrayT& elemmap);
  virtual void ReadGlobalElementMap (const StringT& name, iArrayT& elemmap);
  virtual void ReadGlobalElementSet (const StringT& name, iArrayT& set);
  virtual void ReadConnectivity (const StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (const StringT& name, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet (const StringT& name);
  virtual void ReadNodeSet (const StringT& name, iArrayT& nodes); /* offset nodes, continuous */

  virtual bool AreSideSetsLocal (void) const;
  virtual int  NumSidesInSet (const StringT& setname) const;
  virtual StringT SideSetGroupName (const StringT& setname) const;
  virtual void ReadSideSetLocal (const StringT& setname, iArray2DT& sides) const;
  virtual void ReadSideSetGlobal (const StringT& setname, iArray2DT& sides) const;
  
  virtual void QARecords (ArrayT<StringT>& records);
  virtual int  NumTimeSteps (void) const;
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void) const;
  virtual int  NumElementVariables (void) const;
  virtual int  NumQuadratureVariables (void) const;

  virtual void ReadAllNodeVariable (int step, int varindex, dArrayT& values);
  virtual void ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values);
  virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;  
  
  virtual void NodeVariablesUsed (const StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (const StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (const StringT& name, iArrayT& used);  

  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  virtual void ReadNodeVariables (int step, const StringT& name, dArray2DT& nvalues);
  virtual void ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& nvalues);

  virtual void ReadAllElementVariable (int step, int varindex, dArrayT& values);
  virtual void ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values);
  virtual void ReadAllElementVariables (int step, dArray2DT& evalues);
  virtual void ReadElementVariables (int step, const StringT& name, dArray2DT& evalues);

  virtual void ReadAllQuadratureVariable (int step, int varindex, dArrayT& values);
  virtual void ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values);
  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  virtual void ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues);

 private:
  bool AdvanceStream (istream& in, const char* key) const;
  void ScanGeometryFile (void);
  
  StringT CreateVariableFile (const StringT& old, int inc) const;
  int  ComponentIndex (int varindex, ArrayT<StringT>& labels) const;
  void VarPrelims_Geo (int step, const StringT& name, int& group_id, int& currentinc);
  void VarPrelims_Case (AutoArrayT<bool>& vector, AutoArrayT<StringT>& labels, bool node);
  void ReadOneVariableData (int component, const StringT& label, int group_id, dArrayT& values, int currentinc, bool nodal) const;
  void ReadVariableData (ArrayT<bool>& vector, ArrayT<StringT>& labels, int group_id, dArray2DT& values, int currentinc, bool nodal) const;
  void VariableUsed (const StringT& name, iArrayT& used, ArrayT<StringT>& labels, ArrayT<bool>& vector, bool nodal) const;  

  enum PartDimensionsT { kNumNodes=0, kNumElements, kPartID, kPartDims };

 private:
  EnSightT fData;
  StringT fGeometryFile;
  StringT fCaseFile;
  iArray2DT fPartDimensions; // num_nodes, num_elems, partID
  int fStartIncrement;
  int fIncrement;
};

inline void EnSightInputT::Close (void) { }
inline void EnSightInputT::SideSetNames (ArrayT<StringT>& sidenames) const
{ sidenames.Free (); }
inline void EnSightInputT::NodeSetNames (ArrayT<StringT>& nodenames) const
{ nodenames.Free (); }
inline int EnSightInputT::NumElementQuadPoints (const StringT& name)
{
#pragma unused (name)
  return (0);
}
inline int EnSightInputT::NumSideSets (void) const { return 0; }
inline int EnSightInputT::NumNodeSets (void) const { return 0; }
inline int EnSightInputT::NumDimensions (void) const { return 3; }
inline int EnSightInputT::NumNodesInSet (const StringT& name) 
{ 
#pragma unused(name)
	return 0; 
}
inline void EnSightInputT::ReadNodeSet (const StringT& name, iArrayT& nodes)
{
#pragma unused (name)
  nodes.Free ();
}
inline bool EnSightInputT::AreSideSetsLocal (void) const { return true; }
inline int  EnSightInputT::NumSidesInSet (const StringT& setname) const
{
#pragma unused (setname)
  return 0;
}
inline StringT EnSightInputT::SideSetGroupName (const StringT& setname) const
{
#pragma unused (setname)
  StringT name ("");
  return name; 
}
inline void EnSightInputT::ReadSideSetLocal (const StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free ();
}
inline void EnSightInputT::ReadSideSetGlobal (const StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free ();
}
inline int EnSightInputT::NumQuadratureVariables (void) const { return 0; }
inline void EnSightInputT::QuadratureVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void EnSightInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free (); }
inline void EnSightInputT::ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (nsetname)
  nvalues.Free();
}
inline void EnSightInputT::ReadAllQuadratureVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
}
inline void EnSightInputT::ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
}
inline void EnSightInputT::ReadAllQuadratureVariables (int step, dArray2DT& qvalues)
{
#pragma unused (step)
  qvalues.Free();
}
inline void EnSightInputT::ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues)
{
#pragma unused (step)
#pragma unused (name)
  qvalues.Free();
}

} // namespace Tahoe 
#endif
