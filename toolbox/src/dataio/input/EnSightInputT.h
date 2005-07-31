/* $Id: EnSightInputT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1998)                                          */

#ifndef _ENSIGHTINPUT_T_H_
#define _ENSIGHTINPUT_T_H_

#include "InputBaseT.h"

/* direct members */
#include "EnSightT.h"
#include "StringT.h"
#include "iArrayT.h"
#include "iArray2DT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
template <class TYPE> class ArrayT;
class dArray2DT;
class ifstreamT;

class EnSightInputT : public InputBaseT
{
public:
EnSightInputT (ostream& out, bool binary, const char* filename);

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

// EnSight does not stores side or node sets
inline int EnSightInputT::NumSideSets (void) const { return 0; }
inline int EnSightInputT::NumNodeSets (void) const { return 0; }
inline void EnSightInputT::SideSetNumbers (iArrayT& sidenums) const
{ sidenums.Allocate (0);}
inline void EnSightInputT::NodeSetNumbers (iArrayT& nodenums) const
{ nodenums.Allocate (0); }
inline void EnSightInputT::ReadNodeSet (int set_num, iArrayT& nodes) const
{
#pragma unused (set_num)
#pragma unused (nodes)
}
inline void EnSightInputT::ReadSideSet (int set_num, iArray2DT& sides) const
{
#pragma unused (set_num)
#pragma unused (sides)
}
inline void EnSightInputT::ReadSideSetGlobal (int set_num, iArray2DT& sides) const
{
#pragma unused (set_num)
#pragma unused (sides)
}

#endif
