/* $Id: TahoeInputT.h,v 1.7 2001-12-16 23:53:46 paklein Exp $ */
/* created: sawimme July 2001 */

#ifndef _TAHOEINPUT_T_H_
#define _TAHOEINPUT_T_H_

#include "InputBaseT.h"
#include "ModelFileT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

class TahoeInputT : public InputBaseT
{
 public:
  TahoeInputT (ostream& out);
  
  virtual void Open (const StringT& filename);
  virtual void Close (void);

  virtual void ElementGroupNames (ArrayT<StringT>& groupnames) const;
  virtual void SideSetNames (ArrayT<StringT>& sidenames) const;
  virtual void NodeSetNames (ArrayT<StringT>& nodenames) const;

  virtual int NumElementGroups (void) const;
  virtual int NumSideSets (void) const;
  virtual int NumNodeSets (void) const;

  virtual int NumNodes (void) const;
  virtual int NumDimensions (void) const;
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  virtual int NumGlobalElements (void) const;
  virtual int NumElements (StringT& name);
  virtual int NumElementNodes (StringT& name);
  virtual int  NumElementQuadPoints (StringT& name);
  virtual void ReadAllElementMap (iArrayT& elemmap);
  virtual void ReadGlobalElementMap (StringT& name, iArrayT& elemmap);
  virtual void ReadGlobalElementSet (StringT& name, iArrayT& set);
  virtual void ReadConnectivity (StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (StringT& name, GeometryT::CodeT& code);

  virtual int  NumNodesInSet (StringT& name);
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes);

  virtual bool AreSideSetsLocal (void) const;
  virtual int  NumSidesInSet (StringT& name) const;
  virtual StringT SideSetGroupName (StringT& name) const;
  virtual void ReadSideSetLocal (StringT& name, iArray2DT& sides) const;
  virtual void ReadSideSetGlobal (StringT& name, iArray2DT& sides) const;

  virtual void QARecords (ArrayT<StringT>& records);
  virtual int  NumTimeSteps (void) const;
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void) const;
  virtual int  NumElementVariables (void) const;
  virtual int  NumQuadratureVariables (void) const;
  virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;  
  virtual void NodeVariablesUsed (StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (StringT& name, iArrayT& used);  
  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  virtual void ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues);
  virtual void ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues);
  virtual void ReadAllElementVariables (int step, dArray2DT& evalues);
  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);
  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues);


 private:
  void SetCode (int numelemnodes, int dof, GeometryT::CodeT& code) const;

 private:
  ModelFileT fModel;
};

inline bool TahoeInputT::AreSideSetsLocal (void) const
{ return true; }

inline int TahoeInputT::NumElementQuadPoints (StringT& name)
{
#pragma unused (name)
  return (0);
}
inline void TahoeInputT::QARecords (ArrayT<StringT>& records)
{
//TEMP
#pragma unused(records)
}

inline int  TahoeInputT::NumTimeSteps (void) const
{ return 0; }
inline void TahoeInputT::ReadTimeSteps (dArrayT& steps)
{ steps.Free (); }
inline int  TahoeInputT::NumNodeVariables (void) const
{ return 0; }
inline int  TahoeInputT::NumElementVariables (void) const
{ return 0; }
inline int  TahoeInputT::NumQuadratureVariables (void) const
{ return 0; }
inline void TahoeInputT::NodeVariablesUsed (StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void TahoeInputT::ElementVariablesUsed (StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void TahoeInputT::QuadratureVariablesUsed (StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void TahoeInputT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{ nlabels.Free (); }
inline void TahoeInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{ elabels.Free (); }
inline void TahoeInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const 
{ qlabels.Free (); }
inline void TahoeInputT::ReadAllNodeVariables (int step, dArray2DT& nvalues)
{
#pragma unused (step)
  nvalues.Free (); 
}
inline void TahoeInputT::ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (name)
  nvalues.Free (); 
}
inline void TahoeInputT::ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (nsetname)
  nvalues.Free (); 
}
inline void TahoeInputT::ReadAllElementVariables (int step, dArray2DT& evalues)
{
#pragma unused (step)
  evalues.Free (); 
}
inline void TahoeInputT::ReadElementVariables (int step, StringT& name, dArray2DT& evalues)
{
#pragma unused (step)
#pragma unused (name)
  evalues.Free (); 
}
inline void TahoeInputT::ReadAllQuadratureVariables (int step, dArray2DT& qvalues)
{
#pragma unused (step)
  qvalues.Free (); 
}
inline void TahoeInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues)
{
#pragma unused (step)
#pragma unused (name)
  qvalues.Free (); 
}
#endif
