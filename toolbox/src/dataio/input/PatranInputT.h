/* $Id: PatranInputT.h,v 1.2 2001-08-07 23:11:54 paklein Exp $ */
/* created: sawimme July 2001 */

#ifndef _PATRANINPUT_T_H_
#define _PATRANINPUT_T_H_

#include "InputBaseT.h"
#include "PatranT.h"

class PatranInputT : public InputBaseT
{
 public:
  PatranInputT (ostream& out);
  
  virtual void Open (const StringT& file);
  virtual void Close (void);

  virtual void ElementGroupNames (ArrayT<StringT>& groupnames);
  virtual void SideSetNames (ArrayT<StringT>& sidenames);
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
  virtual int  NumSidesInSet (StringT& name);
  virtual int  SideSetGroupIndex (StringT& name);
  virtual void ReadSideSetLocal (StringT& name, iArray2DT& sides);
  virtual void ReadSideSetGlobal (StringT& name, iArray2DT& sides);

 private:
  void SetCode (int namedtype, GeometryT::CodeT& code) const;

 private:
  PatranT fPatran;
};

inline void PatranInputT::SideSetNames (ArrayT<StringT>& sidenames)
{ 
#pragma unused (sidenames) 
}
inline int  PatranInputT::NumSideSets (void)
{ return 0; }
inline int PatranInputT::NumNodes (void)
{ return fPatran.NumNodes(); }
inline int PatranInputT::NumDimensions (void)
{ return fPatran.NumDimensions(); }
inline int PatranInputT::NumGlobalElements (void)
{ return fPatran.NumElements (); }
inline bool PatranInputT::AreSideSetsLocal (void)
{ return false; }




#endif
