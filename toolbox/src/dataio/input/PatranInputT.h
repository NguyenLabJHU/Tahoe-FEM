/* $Id: PatranInputT.h,v 1.12 2002-07-25 19:47:29 sawimme Exp $ */
/* created: sawimme July 2001 */

#ifndef _PATRANINPUT_T_H_
#define _PATRANINPUT_T_H_

#include "InputBaseT.h"
#include "PatranT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArrayT.h"


namespace Tahoe {

class PatranInputT : public InputBaseT
{
 public:
  PatranInputT (ostream& out);
  
  virtual bool Open (const StringT& file);
  virtual void Close (void);

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
  virtual int  NumSidesInSet (const StringT& name) const;
  virtual StringT SideSetGroupName (const StringT& name) const;
  virtual void ReadSideSetLocal (const StringT& name, iArray2DT& sides) const;
  virtual void ReadSideSetGlobal (const StringT& name, iArray2DT& sides) const;

  virtual void QARecords (ArrayT<StringT>& records);
  virtual int  NumTimeSteps (void) const;
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void) const;
  virtual int  NumElementVariables (void) const;
  virtual int  NumQuadratureVariables (void) const;

  virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;  

  virtual void NodeVariablesUsed (const StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (const StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (const StringT& name, iArrayT& used);  

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
  void SetCode (PatranT::NamedTypes namedtype, GeometryT::CodeT& code) const;

 private:
  PatranT fPatran;
};

inline void PatranInputT::SideSetNames (ArrayT<StringT>& sidenames) const
{ 
#pragma unused (sidenames) 
}
inline int  PatranInputT::NumSideSets (void) const
{ return 0; }
inline int PatranInputT::NumNodes (void) const
{ return fPatran.NumNodes(); }
inline int PatranInputT::NumDimensions (void) const
{ return fPatran.NumDimensions(); }
inline int PatranInputT::NumElementQuadPoints (const StringT& name)
{
#pragma unused (name)
  return (0);
}
inline int PatranInputT::NumGlobalElements (void) const
{ return fPatran.NumElements (); }
inline bool PatranInputT::AreSideSetsLocal (void) const
{ return false; }
inline void PatranInputT::QARecords (ArrayT<StringT>& records) 
{ fPatran.VersionNotes (records); }
inline int  PatranInputT::NumTimeSteps (void) const
{ return 0; }
inline void PatranInputT::ReadTimeSteps (dArrayT& steps)
{ steps.Free (); }
inline int  PatranInputT::NumNodeVariables (void) const
{ return 0; }
inline int  PatranInputT::NumElementVariables (void) const
{ return 0; }
inline int  PatranInputT::NumQuadratureVariables (void) const
{ return 0; }
inline void PatranInputT::NodeVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void PatranInputT::ElementVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void PatranInputT::QuadratureVariablesUsed (const StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void PatranInputT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{ nlabels.Free (); }
inline void PatranInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{ elabels.Free (); }
inline void PatranInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free (); } 
inline void PatranInputT::ReadAllNodeVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
}
inline void PatranInputT::ReadNodeVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
}
inline void PatranInputT::ReadAllNodeVariables (int step, dArray2DT& nvalues)
{
#pragma unused (step)
  nvalues.Free ();
}
inline void PatranInputT::ReadNodeVariables (int step, const StringT& name, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (name)
  nvalues.Free ();
}
inline void PatranInputT::ReadNodeSetVariables (int step, const StringT& nsetname, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (nsetname)
  nvalues.Free ();
}
inline void PatranInputT::ReadAllElementVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
}
inline void PatranInputT::ReadElementVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
}
inline void PatranInputT::ReadAllElementVariables (int step, dArray2DT& evalues)
{
#pragma unused (step)
  evalues.Free ();
}
inline void PatranInputT::ReadElementVariables (int step, const StringT& name, dArray2DT& evalues)
{
#pragma unused (step)
#pragma unused (name)
  evalues.Free ();
}
inline void PatranInputT::ReadAllQuadratureVariable (int step, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (varindex)
  values.Free();
}
inline void PatranInputT::ReadQuadratureVariable (int step, const StringT& name, int varindex, dArrayT& values)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (varindex)
  values.Free();
}
inline void PatranInputT::ReadAllQuadratureVariables (int step, dArray2DT& qvalues)
{
#pragma unused (step)
  qvalues.Free ();
}
inline void PatranInputT::ReadQuadratureVariables (int step, const StringT& name, dArray2DT& qvalues)
{
#pragma unused (step)
#pragma unused (name)
  qvalues.Free ();
}


} // namespace Tahoe 
#endif
