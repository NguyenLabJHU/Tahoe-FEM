#ifndef _INPUTFEASCII_T_H_
#define _INPUTFEASCII_T_H_

#include "InputBaseT.h"
#include "StringT.h"
#include "iAutoArrayT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"
class ifstreamT;
class dArrayT;

class InputFEASCIIT : public InputBaseT
{
public:
  InputFEASCIIT (ostream& out);

  virtual void Open (const StringT& filename);
  virtual void Close (void);

  /* return names, Array must be preallocated */
  virtual void ElementGroupNames (ArrayT<StringT>& groupnames) const;
  virtual void SideSetNames (ArrayT<StringT>& sidenames) const;
  virtual void NodeSetNames (ArrayT<StringT>& nodenames) const;

  /* return dimenesions */
  virtual int  NumElementGroups (void) const;
  virtual int  NumSideSets (void) const;
  virtual int  NumNodeSets (void) const;

  /* NODES */
  virtual int  NumNodes (void) const;
  virtual int  NumDimensions (void) const;
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  /* ELEMENTS */
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
  
  /* record[0] = progname, record[1] = version, record[2] = date, record[3] = time */
  virtual void QARecords (ArrayT<StringT>& records);

  virtual int  NumTimeSteps (void) const;
  virtual void ReadTimeSteps (dArrayT& steps);

  /* for all nodes or elements */
  virtual int  NumNodeVariables (void) const;
  virtual int  NumElementVariables (void) const;
  virtual int  NumQuadratureVariables (void) const;

  virtual void NodeVariablesUsed (StringT& name, iArrayT& used);
  virtual void ElementVariablesUsed (StringT& name, iArrayT& used);
  virtual void QuadratureVariablesUsed (StringT& name, iArrayT& used);  

  virtual void ReadNodeLabels (ArrayT<StringT>& labels) const;
  virtual void ReadElementLabels (ArrayT<StringT>& elabels) const;
  virtual void ReadQuadratureLabels (ArrayT<StringT>& qlabels) const;

  /* step starts at zero and increases by one */
  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  virtual void ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues);
  virtual void ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues);

  virtual void ReadAllElementVariables (int step, dArray2DT& evalues);
  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);

  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues);
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues);

 private:
  bool OpenFile (ifstreamT& in, const char *ext) const;
  bool ScanGeometryFile (ifstreamT& in);
  bool ScanResultsFile (ifstreamT& in);
  bool AdvanceToBlock (ifstreamT& in, const StringT& name, const char *tname) const;
  void DataBlock (ifstreamT& in, iArrayT& used, iArrayT& ids, dArray2DT& vals, bool nodal) const;

 private:
  StringT fFileRoot;

  iAutoArrayT fBlockID;
  int fNumNodes;
  int fNumElements;
  int fNumDOF;

  AutoArrayT<double> fTimeSteps;
  AutoArrayT<StringT> fNodeVariable;
  AutoArrayT<StringT> fElementVariable;
};

inline void InputFEASCIIT::SideSetNames (ArrayT<StringT>& sidenames) const
{ sidenames.Free(); }
inline void InputFEASCIIT::NodeSetNames (ArrayT<StringT>& nodenames) const
{ nodenames.Free(); }
inline int InputFEASCIIT::NumElementGroups (void) const
{ return fBlockID.Length(); }
inline int InputFEASCIIT::NumSideSets (void) const
{ return 0; }
inline int InputFEASCIIT::NumNodeSets (void) const
{ return 0; }
inline int InputFEASCIIT::NumNodes (void) const
{ return fNumNodes; }
inline int InputFEASCIIT::NumDimensions (void) const
{ return fNumDOF; }
inline int InputFEASCIIT::NumElementQuadPoints (StringT& name)
{ return 0; }
inline int InputFEASCIIT::NumGlobalElements (void) const
{ return fNumElements; }
inline int InputFEASCIIT::NumNodesInSet (StringT& name)
{
#pragma unused (name)
  return 0;
}
inline void InputFEASCIIT::ReadNodeSet (StringT& name, iArrayT& nodes)
{
#pragma unused (name)
  nodes.Free();
}
inline bool InputFEASCIIT::AreSideSetsLocal (void) const
{ return true; }
inline int InputFEASCIIT::NumSidesInSet (StringT& setname) const
{ 
#pragma unused (setname)
  return 0; 
}
inline StringT InputFEASCIIT::SideSetGroupName (StringT& setname) const
{ 
#pragma unused (setname)
  StringT s = "";
  return s;
}
inline void InputFEASCIIT::ReadSideSetLocal (StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free ();
}
inline void InputFEASCIIT::ReadSideSetGlobal (StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free ();
}
inline void InputFEASCIIT::QARecords (ArrayT<StringT>& records)
{
#pragma unused (records)
}
inline int InputFEASCIIT::NumTimeSteps (void) const
{ return fTimeSteps.Length(); }
inline int InputFEASCIIT::NumNodeVariables (void) const
{ return fNodeVariable.Length(); }
inline int InputFEASCIIT::NumElementVariables (void) const
{ return fElementVariable.Length(); }
inline int InputFEASCIIT::NumQuadratureVariables (void) const
{ return 0; }
inline void InputFEASCIIT::QuadratureVariablesUsed (StringT& name, iArrayT& used)
{
#pragma unused (name)
  used = 0;
}
inline void InputFEASCIIT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free(); }
inline void InputFEASCIIT::ReadNodeSetVariables (int step, StringT& name, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (name)
  nvalues.Free();
}
inline void InputFEASCIIT::ReadAllQuadratureVariables (int step, dArray2DT& vals)
{
#pragma unused (step)
  vals.Free();
}
inline void InputFEASCIIT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& vals)
{
#pragma unused (step)
#pragma unused (name)
  vals.Free();
}

#endif
