/* $Id: ExodusInputT.h,v 1.7 2001-12-16 23:53:45 paklein Exp $ */
/* created: sawimme (05/18/1998) */

#ifndef _EXODUSINPUT_T_H_
#define _EXODUSINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "ExodusT.h"
#include "dArray2DT.h"
#include "iArrayT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class iArray2DT;

class ExodusInputT : public InputBaseT
{
public:
  ExodusInputT (ostream& out);

  virtual void Open (const StringT& file);
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
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  virtual int  NumGlobalElements (void) const;
  virtual int  NumElements (StringT& name);
  virtual int  NumElementNodes (StringT& name);
  virtual int  NumElementQuadPoints (StringT& name);
  virtual void ReadAllElementMap (iArrayT& elemmap);
  virtual void ReadGlobalElementMap (StringT& name, iArrayT& elemmap);
  virtual void ReadGlobalElementSet (StringT& name, iArrayT& set);
  virtual void ReadConnectivity (StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet (StringT& name);
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes);

  virtual bool AreSideSetsLocal (void) const;
  virtual int  NumSidesInSet (StringT& setname) const;
  virtual StringT SideSetGroupName (StringT& setname) const;
  virtual void ReadSideSetLocal (StringT& setname, iArray2DT& sides) const;
  virtual void ReadSideSetGlobal (StringT& setname, iArray2DT& sides) const;

  virtual void QARecords (ArrayT<StringT>& records);

  virtual int  NumTimeSteps (void) const;
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void) const;
  virtual int  NumElementVariables (void) const;
  virtual int  NumQuadratureVariables (void) const;

  virtual void NodeVariablesUsed (StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (StringT& name, iArrayT& used);  

  virtual void ReadNodeLabels (ArrayT<StringT>& labels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;

  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  virtual void ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues);
  virtual void ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues);

  virtual void ReadAllElementVariables (int step, dArray2DT& evalues);
  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);

  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues);

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

inline int ExodusInputT::NumElementQuadPoints (StringT& name)
{
#pragma unused (name)
  return (0);
}
inline int ExodusInputT::NumNodesInSet (StringT& name)
{ 
  int setnum = atoi (name.Pointer());
  return fData.NumNodesInSet (setnum); 
}

inline bool ExodusInputT::AreSideSetsLocal (void) const
{ return true; }

inline int ExodusInputT::NumSidesInSet (StringT& name) const
{  
  int setnum = atoi (name.Pointer());
  return fData.NumSidesInSet (setnum); 
}

inline void ExodusInputT::QARecords (ArrayT<StringT>& records)
{ fData.ReadQA (records); }

inline int ExodusInputT::NumTimeSteps (void) const
{ return fData.NumTimeSteps (); }

inline int ExodusInputT::NumNodeVariables (void) const
{ return fData.NumNodeVariables (); }

inline int ExodusInputT::NumElementVariables (void) const
{ return fData.NumElementVariables (); }

inline int ExodusInputT::NumQuadratureVariables (void) const
{ return 0; }

inline void ExodusInputT::QuadratureVariablesUsed (StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}

inline void ExodusInputT::ReadNodeLabels (ArrayT<StringT>& labels) const
{ fData.ReadNodeLabels (labels); }

inline void ExodusInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{ fData.ReadElementLabels (elabels); }

inline void ExodusInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free (); }

inline void ExodusInputT::ReadAllQuadratureVariables (int step, dArray2DT& vals)
{
#pragma unused (step)
  vals.Free ();
}

inline void ExodusInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& vals)
{
#pragma unused (step)
#pragma unused (name)
  vals.Free ();
}

#endif
