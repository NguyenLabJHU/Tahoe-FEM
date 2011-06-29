/* $Id: ExodusInputT.h,v 1.14 2002-07-05 22:26:26 paklein Exp $ */
/* created: sawimme (05/18/1998) */

#ifndef _EXODUSINPUT_T_H_
#define _EXODUSINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArrayT.h"
#include "dArrayT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class iArray2DT;

class ExodusInputT : public InputBaseT
{
public:
  ExodusInputT (ostream& out);

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
  virtual void ReadNodeSet (const StringT& name, iArrayT& nodes);

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

  virtual void NodeVariablesUsed (const StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (const StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (const StringT& name, iArrayT& used);  

  virtual void ReadNodeLabels (ArrayT<StringT>& labels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;

  virtual void ReadAllNodeVariable (int step, int varindex, dArrayT& values);
  virtual void ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values);
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
  void NodesUsed(const nArrayT<int>& connects, iArrayT& nodesused) const;
  
 private:
  ExodusT fData;
};

inline void ExodusInputT::Close (void)
{ fData.Close (); }

inline int ExodusInputT::NumElementGroups (void) const
{ return fData.NumElementBlocks (); }

inline int ExodusInputT::NumSideSets (void) const
{ return fData.NumSideSets (); }

inline int ExodusInputT::NumNodeSets (void) const
{ return fData.NumNodeSets (); }

inline int ExodusInputT::NumNodes (void) const
{ return fData.NumNodes(); }

inline int ExodusInputT::NumDimensions (void) const
{ return fData.NumDimensions (); }

inline int ExodusInputT::NumElementQuadPoints (const StringT& name)
{
#pragma unused (name)
  return (0);
}
inline int ExodusInputT::NumNodesInSet (const StringT& name)
{ 
  int setnum = atoi (name.Pointer());
  return fData.NumNodesInSet (setnum); 
}

inline bool ExodusInputT::AreSideSetsLocal (void) const
{ return true; }

inline int ExodusInputT::NumSidesInSet (const StringT& name) const
{  
  int setnum = atoi (name.Pointer());
  return fData.NumSidesInSet (setnum); 
}

inline void ExodusInputT::QARecords (ArrayT<StringT>& records)
{ fData.ReadQA (records); }

inline int ExodusInputT::NumTimeSteps (void) const
{ return fData.NumTimeSteps (); }

inline int ExodusInputT::NumNodeVariables (void) const
{ return fData.NumVariables (ExodusT::kNode); }

inline int ExodusInputT::NumElementVariables (void) const
{ return fData.NumVariables (ExodusT::kElement); }

inline int ExodusInputT::NumQuadratureVariables (void) const
{ return 0; }

inline void ExodusInputT::QuadratureVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}

inline void ExodusInputT::ReadNodeLabels (ArrayT<StringT>& labels) const
{ fData.ReadLabels (labels, ExodusT::kNode); }

inline void ExodusInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{ fData.ReadLabels (elabels, ExodusT::kElement); }

inline void ExodusInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free (); }

inline void ExodusInputT::ReadAllQuadratureVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
}
inline void ExodusInputT::ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
}
inline void ExodusInputT::ReadAllQuadratureVariables (int step, dArray2DT& vals)
{
#pragma unused (step)
  vals.Free ();
}

inline void ExodusInputT::ReadQuadratureVariables (int step, const StringT& name, dArray2DT& vals)
{
#pragma unused (step)
#pragma unused (name)
  vals.Free ();
}

} // namespace Tahoe 
#endif
