/* $Id: ExodusInputT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
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
ExodusInputT (ostream& out, const char* filename);

/* virtual with InputManager base class */
virtual int  NumElementGroups (void) const;
virtual int  NumSideSets (void) const;
virtual int  NumNodeSets (void) const;
virtual void GroupNumbers (iArrayT& groupnums) const;
virtual void SideSetNumbers (iArrayT& sidenums) const;
virtual void NodeSetNumbers (iArrayT& nodenums) const;

virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap);
virtual void ReadConnectivity (int group, GeometryT::CodeT& geocode, iArray2DT& connects, iArrayT& elementmap);
virtual void ReadNodeSet (int set_num, iArrayT& nodes) const;
virtual void ReadSideSet (int set_num, iArray2DT& sides) const;
virtual void ReadSideSetGlobal (int set_num, iArray2DT& sides) const;
virtual void Close (void);
virtual void QARecords (ArrayT<StringT>& records) const;
virtual void ReadTimeSteps (dArrayT& steps);
virtual void ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id);
virtual void ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues);

private:
void NodesUsed(const iArray2DT& connects, iArrayT& nodesused) const;

private:
ExodusT fData;
};

inline int ExodusInputT::NumElementGroups (void) const
{ return fData.NumElementBlocks (); }

inline int ExodusInputT::NumSideSets (void) const
{ return fData.NumSideSets (); }

inline int ExodusInputT::NumNodeSets (void) const
{ return fData.NumNodeSets (); }

#endif
