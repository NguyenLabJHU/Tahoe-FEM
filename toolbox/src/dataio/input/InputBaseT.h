/* $Id: InputBaseT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (08/12/1999)                                          */

#ifndef _INPUTBASE_T_H_
#define _INPUTBASE_T_H_

#include "IOBaseT.h"

#include "GeometryT.h"

/* foward declaration */
#include "ios_fwd_decl.h"
class iArrayT;
class iArray2DT;
class dArrayT;
class dArray2DT;
class StringT;
template <class TYPE> class ArrayT;

class InputBaseT : public IOBaseT
{
public:
InputBaseT (ostream& out);

virtual ~InputBaseT (void);

/* virtual with derived classes */
virtual int  NumElementGroups (void) const = 0;
virtual int  NumSideSets (void) const = 0;
virtual int  NumNodeSets (void) const = 0;
virtual void GroupNumbers (iArrayT& groupnums) const = 0;
virtual void SideSetNumbers (iArrayT& sidenums) const = 0;
virtual void NodeSetNumbers (iArrayT& nodenums) const = 0;

virtual void ReadCoordinates (dArray2DT& coords, iArrayT& nodemap) = 0;
virtual void ReadConnectivity (int group, GeometryT::CodeT& geocode, iArray2DT& connects,
	iArrayT& elementmap) = 0;
virtual void ReadNodeSet (int set_num, iArrayT& nodes) const = 0;
virtual void ReadSideSet (int set_num, iArray2DT& sides) const = 0;
virtual void ReadSideSetGlobal (int set_num, iArray2DT& sides) const = 0;

virtual void Close (void) = 0;
virtual void QARecords (ArrayT<StringT>& records) const = 0;

virtual void ReadTimeSteps (dArrayT& steps) = 0;
virtual void ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id) = 0;
virtual void ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues) = 0;
};

#endif
