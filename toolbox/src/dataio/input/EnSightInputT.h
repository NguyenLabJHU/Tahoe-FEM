/* $Id: EnSightInputT.h,v 1.4.2.1 2001-10-15 19:01:27 sawimme Exp $ */
/* created: sawimme (05/18/1998)                                          */

#ifndef _ENSIGHTINPUT_T_H_
#define _ENSIGHTINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "EnSightT.h"
#include "StringT.h"
#include "iArray2DT.h"
#include "iArrayT.h"
#include "dArray2DT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
template <class TYPE> class ArrayT;

class iArray2DT;

class EnSightInputT : public InputBaseT
{
public:
  EnSightInputT (ostream& out, bool binary);

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
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes); /* offset nodes, continuous */

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

  virtual void ReadNodeLabels (ArrayT<StringT>& nlabels) const;
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
  bool AdvanceStream (istream& in, const char* key) const;
  void ScanGeometryFile (void);
  
  StringT CreateVariableFile (const StringT& old, int inc) const;
  void ReadVariableData (ArrayT<bool>& vector, ArrayT<StringT>& labels, int group_id, dArray2DT& values, int currentinc, bool nodal) const;
  
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
inline int EnSightInputT::NumElementQuadPoints (StringT& name)
{
#pragma unused (name)
  return (0);
}
inline int EnSightInputT::NumSideSets (void) const { return 0; }
inline int EnSightInputT::NumNodeSets (void) const { return 0; }
inline int EnSightInputT::NumDimensions (void) const { return 3; }
inline int EnSightInputT::NumNodesInSet (StringT& name) { return 0; }
inline void EnSightInputT::ReadNodeSet (StringT& name, iArrayT& nodes)
{
#pragma unused (name)
  nodes.Free ();
}
inline bool EnSightInputT::AreSideSetsLocal (void) const { return true; }
inline int  EnSightInputT::NumSidesInSet (StringT& setname) const
{
#pragma unused (setname)
  return 0;
}
inline StringT EnSightInputT::SideSetGroupName (StringT& setname) const
{
#pragma unused (setname)
  StringT name ("");
  return name; 
}
inline void EnSightInputT::ReadSideSetLocal (StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free ();
}
inline void EnSightInputT::ReadSideSetGlobal (StringT& setname, iArray2DT& sides) const
{
#pragma unused (setname)
  sides.Free ();
}
inline int EnSightInputT::NumQuadratureVariables (void) const { return 0; }
inline void EnSightInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free (); }
inline void EnSightInputT::ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (nsetname)
  nvalues.Free();
}
inline void EnSightInputT::ReadAllQuadratureVariables (int step, dArray2DT& qvalues)
{
#pragma unused (step)
  qvalues.Free();
}
inline void EnSightInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues)
{
#pragma unused (step)
#pragma unused (name)
  qvalues.Free();
}

#endif
