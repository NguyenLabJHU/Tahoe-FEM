/* $Id: PatranInputT.h,v 1.3 2001-09-04 14:46:38 sawimme Exp $ */
/* created: sawimme July 2001 */

#ifndef _PATRANINPUT_T_H_
#define _PATRANINPUT_T_H_

#include "InputBaseT.h"
#include "PatranT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

class PatranInputT : public InputBaseT
{
 public:
  PatranInputT (ostream& out);
  
  virtual void Open (const StringT& file);
  virtual void Close (void);

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
  virtual void ReadAllElementMap (iArrayT& elemmap);
  virtual void ReadGlobalElementMap (StringT& name, iArrayT& elemmap);
  virtual void ReadGlobalElementSet (StringT& name, iArrayT& set);
  virtual void ReadConnectivity (StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet (StringT& name);
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes);

  virtual bool AreSideSetsLocal (void) const;
  virtual int  NumSidesInSet (StringT& name) const;
  virtual int  SideSetGroupIndex (StringT& name) const;
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

  virtual void ReadAllNodeVariables (int step, dArray2DT& nvalues);
  virtual void ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues);
  virtual void ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues); 

  virtual void ReadAllElementVariables (int step, dArray2DT& evalues); 
  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);

  virtual void ReadAllQuadratureVariables (int step, dArray2DT& qvalues); 
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues); 

 private:
  void SetCode (int namedtype, GeometryT::CodeT& code) const;

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
inline void PatranInputT::ReadNodeLabels (ArrayT<StringT>& nlabels) const
{ nlabels.Free (); }
inline void PatranInputT::ReadElementLabels (ArrayT<StringT>& elabels) const
{ elabels.Free (); }
inline void PatranInputT::ReadQuadratureLabels (ArrayT<StringT>& qlabels) const
{ qlabels.Free (); } 
inline void PatranInputT::ReadAllNodeVariables (int step, dArray2DT& nvalues)
{
#pragma unused (step)
  nvalues.Free ();
}
inline void PatranInputT::ReadNodeVariables (int step, StringT& name, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (name)
  nvalues.Free ();
}
inline void PatranInputT::ReadNodeSetVariables (int step, StringT& nsetname, dArray2DT& nvalues)
{
#pragma unused (step)
#pragma unused (nsetname)
  nvalues.Free ();
}
inline void PatranInputT::ReadAllElementVariables (int step, dArray2DT& evalues)
{
#pragma unused (step)
  evalues.Free ();
}
inline void PatranInputT::ReadElementVariables (int step, StringT& name, dArray2DT& evalues)
{
#pragma unused (step)
#pragma unused (name)
  evalues.Free ();
}
inline void PatranInputT::ReadAllQuadratureVariables (int step, dArray2DT& qvalues)
{
#pragma unused (step)
  qvalues.Free ();
}
inline void PatranInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues)
{
#pragma unused (step)
#pragma unused (name)
  qvalues.Free ();
}


#endif
