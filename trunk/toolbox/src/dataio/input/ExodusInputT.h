/* $Id: ExodusInputT.h,v 1.2 2001-08-03 19:16:43 sawimme Exp $ */
/* created: sawimme (05/18/1998)                                          */

#ifndef _EXODUSINPUT_T_H_
#define _EXODUSINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "ExodusT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
template <class TYPE> class ArrayT;
class dArray2DT;
class iArray2DT;

class ExodusInputT : public InputBaseT
{
public:
  ExodusInputT (ostream& out);

  virtual void Open (const StringT& file);
  virtual void Close (void);

  /* virtual with InputManager base class */
  virtual int  NumElementGroups (void);
  virtual int  NumSideSets (void);
  virtual int  NumNodeSets (void);

  virtual int  NumNodes (void);
  virtual int  NumDimensions (void);
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  virtual bool AreSideSetsLocal (void);

  virtual int  NumGlobalElements (void);
  virtual void ReadAllElementMap (iArrayT& elemmap);

  virtual void QARecords (ArrayT<StringT>& records);

  virtual int  NumTimeSteps (void);
  virtual void ReadTimeSteps (dArrayT& steps);

  virtual int  NumNodeVariables (void);
  virtual int  NumElementVariables (void);

 protected:
  virtual void ElementGroupIDs (iArrayT& groupnums);
  virtual void SideSetIDs (iArrayT& sidenums);
  virtual void NodeSetIDs (iArrayT& nodenums);
  
  virtual int NumElements_ID (int ID);
  virtual int NumElementNodes_ID (int ID);
  virtual void ReadGlobalElementMap_ID (int ID, iArrayT& elemmap);
  virtual void ReadConnectivity_ID (int ID, iArray2DT& connects);
  virtual void ReadGeometryCode_ID (int ID, GeometryT::CodeT& code);

  virtual int  NumNodesInSet_ID (int ID);
  virtual void ReadNodeSet_ID (int ID, iArrayT& nodes);

  virtual int  NumSidesInSet_ID (int ID);
  virtual int  SideSetGroupIndex_ID (int ID);
  virtual void ReadSideSetLocal_ID (int ID, iArray2DT& sides);
  virtual void ReadSideSetGlobal_ID (int ID, iArray2DT& sides);

  virtual void ReadElementLabels_ID (int ID, ArrayT<StringT>& elabels);
  virtual void ReadElementVariables_ID (int step, int ID, dArray2DT& evalues);

 private:
  void NodesUsed(const iArray2DT& connects, iArrayT& nodesused) const;
  
 private:
  ExodusT fData;
};

inline int ExodusInputT::NumDimensions (void)
{ return fData.NumDimensions (); }

inline int ExodusInputT::NumElementGroups (void)
{ return fData.NumElementBlocks (); }

inline int ExodusInputT::NumSideSets (void)
{ return fData.NumSideSets (); }

inline int ExodusInputT::NumNodeSets (void)
{ return fData.NumNodeSets (); }

inline int ExodusInputT::NumNodes (void)
{ return fData.NumNodes(); }

inline void ExodusInputT::Close (void)
{ fData.Close (); }

inline void ExodusInputT::QARecords (ArrayT<StringT>& records)
{ fData.ReadQA (records); }

inline int ExodusInputT::NumTimeSteps (void)
{ return fData.NumTimeSteps (); }

inline bool ExodusInputT::AreSideSetsLocal (void)
{ return true; }

inline int ExodusInputT::NumNodeVariables (void)
{ return fData.NumNodeVariables (); }

inline int ExodusInputT::NumElementVariables (void)
{ return fData.NumElementVariables (); }

inline int ExodusInputT::NumNodesInSet_ID (int setnum)
{ return fData.NumNodesInSet (setnum); }

inline int ExodusInputT::NumSidesInSet_ID (int setnum)
{ return fData.NumSidesInSet (setnum); }

#endif
