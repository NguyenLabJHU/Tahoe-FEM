/* $Id: AbaqusInputT.h,v 1.3 2001-08-08 11:47:12 sawimme Exp $ */
/* created: sawimme (05/18/1998)                                          */

#ifndef _ABAQUSINPUT_T_H_
#define _ABAQUSINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "AbaqusT.h"
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
  AbaqusInputT (ostream& out, bool binary);

  virtual void Open (const StringT& file);
  virtual void Close (void);

  virtual void ElementGroupNames (ArrayT<StringT>& groupnames);
  virtual void NodeSetNames (ArrayT<StringT>& nodenames);

  virtual int  NumElementGroups (void);
  virtual int  NumSideSets (void);
  virtual int  NumNodeSets (void);

  virtual int  NumNodes (void);
  virtual int  NumDimensions (void);
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  virtual int  NumGlobalElements (void);
  virtual int  NumElements (StringT& name);
  virtual int  NumElementNodes (StringT& name);
  virtual void ReadAllElementMap (iArrayT& elemmap);
  virtual void ReadGlobalElementMap (StringT& name, iArrayT& elemmap);
  virtual void ReadConnectivity (StringT& name, iArray2DT& connects);
  virtual void ReadGeometryCode (StringT& name, GeometryT::CodeT& geocode);

  virtual int  NumNodesInSet (StringT& name);
  virtual void ReadNodeSet (StringT& name, iArrayT& nodes);

  virtual bool AreSideSetsLocal (void);

  virtual void QARecords (ArrayT<StringT>& records);

  virtual int  NumTimeSteps (void);
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void);
  virtual int  NumElementVariables (void);
  virtual int  NumQuadratureVariables (void);

  virtual void ReadElementLabels (StringT& name, ArrayT<StringT>& elabels);
  virtual void ReadQuadratureLabels (StringT& name, ArrayT<StringT>& qlabels);  

  virtual void ReadElementVariables (int step, StringT& name, dArray2DT& evalues);
  virtual void ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues);

private:
  bool OpenFile (ifstream& in);
  void Initialize (void);
  void ReadElementSets (void);
  void ReadNodeSets (void);
  void SetLabelName (AbaqusT::VariableKeyT key, int& unknown, StringT& name, char incrementor) const;
  
 private:
  AbaqusT fData;
  StringT fResultFile;
  StringT fVersion, fDate, fTime;
  int fNumElements;
  int fNumNodes;
  int fDimensions;
  
  iArray2DT fElementData; // elem ID, geocode, setname index, numelemnodes, numelems
  ArrayT<StringT> fElementSetNames;
  iArray2DT fCoordinateData; // node ID, setname index
  ArrayT<StringT> fNodeSetNames;
  
  // variable data
  bool fEnergyData;
};

inline int AbaqusInputT::NumElementGroups (void) { return fElementSetNames.Length(); }
inline int AbaqusInputT::NumNodeSets (void) { return fNodeSetNames.Length(); }
inline int AbaqusInputT::NumNodes (void) { return fNumNodes; }
inline int AbaqusInputT::NumDimensions (void) { return fDimensions; }
inline int AbaqusInputT::NumGlobalElements (void) { return fNumElements; }
inline bool AbaqusInputT::AreSideSetsLocal (void) { return false; }
inline int AbaqusInputT::NumQuadratureVariables (void) { return 0; } // add later if needed
inline void AbaqusInputT::ReadQuadratureLabels (StringT& name, ArrayT<StringT>& qlabels)
{
#pragma unused (name)
#pragma unused (qlabels)
}
inline void AbaqusInputT::ReadQuadratureVariables (int step, StringT& name, dArray2DT& qvalues)
{
#pragma unused (step)
#pragma unused (name)
#pragma unused (qvalues)
}

#endif
