/* created: sawimme July 2001 */

#ifndef _TAHOEINPUT_T_H_
#define _TAHOEINPUT_T_H_

#include "InputBaseT.h"
#include "ModelFileT.h"

class TahoeInputT : public InputBaseT
{
 public:
  TahoeInputT (ostream& out);
  
  virtual void Open (const StringT& filename);
  virtual void Close (void);

  virtual int NumElementGroups (void);
  virtual int NumSideSets (void);
  virtual int NumNodeSets (void);

  virtual int NumNodes (void);
  virtual int NumDimensions (void);
  virtual void ReadNodeMap (iArrayT& nodemap);
  virtual void ReadCoordinates (dArray2DT& coords);
  virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);

  virtual bool AreSideSetsLocal (void);

  virtual int NumGlobalElements (void);
  virtual void ReadAllElementMap (iArrayT& elemmap);

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

 private:
  void SetCode (int numelemnodes, int dof, GeometryT::CodeT& code);

 private:
  ModelFileT fModel;
};

inline bool TahoeInputT::AreSideSetsLocal (void)
{ return true; }

#endif
