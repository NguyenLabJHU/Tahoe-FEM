/* $Id: AbaqusInputT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1998)                                          */
/* If SetLabelName does not have a case for the variable being read, a default variable name is used. As new variables are added to AbaqusT::VariableKey, they should also be added to SetLabelName */

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

class AbaqusInputT : public InputBaseT
{
public:
AbaqusInputT (ostream& out, bool binary, const char* filename);

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
int fNumElementSets;
int fNumNodeSets;

iArray2DT fElementData; // elem ID, geocode, set ID, numelemnodes, numelems
iArray2DT fCoordinateData; // node ID, set ID

// variable data
bool fEnergyData;
};

// Abaqus does not stores side sets
inline int AbaqusInputT::NumSideSets (void) const { return 0; }
inline void AbaqusInputT::SideSetNumbers (iArrayT& sidenums) const
{ sidenums.Allocate (0);}
inline void AbaqusInputT::ReadSideSet (int set_num, iArray2DT& sides) const
{
#pragma unused (set_num)
#pragma unused (sides)
}
inline void AbaqusInputT::ReadSideSetGlobal (int set_num, iArray2DT& sides) const
{
#pragma unused (set_num)
#pragma unused (sides)
}

#endif
