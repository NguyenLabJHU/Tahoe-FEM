/* $Id: EnSightOutputT.h,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: sawimme (05/18/1999)                                          */

#ifndef _ENSIGHTOUTPUT_T_H_
#define _ENSIGHTOUTPUT_T_H_

#include "OutputBaseT.h"
#include "AutoArrayT.h"
#include "StringT.h"
#include "EnSightT.h"

class EnSightOutputT : public OutputBaseT
{
public:
EnSightOutputT (ostream& out, const ArrayT<StringT>& out_strings,
	int numdigs, bool binary);

void WriteGeometry (void);
void WriteOutput (double time, int ID, const dArray2DT& n_values,
	const dArray2DT& e_values);

private:

enum FileNameTypeT { kWildFile = -9, kNoIncFile = -1 };

StringT OpenGeometryFile (EnSightT& ens, ofstream& geo) const;
StringT CreateFileName (const StringT& label, int increment) const;

void WritePart (ostream& geo, EnSightT& ens, int ID) const;
void WriteCoordinates (ostream& geo, EnSightT& ens, const iArrayT& nodes_used) const;
void WriteConnectivity (ostream& geo, EnSightT& ens, const iArrayT& nodes_used, int index) const;

void WriteVariable (EnSightT& ens, bool nodal, int ID, const dArray2DT& values,
	const ArrayT<StringT>& labels, AutoArrayT<StringT>& names,
	AutoArrayT<StringT>& files, AutoArrayT<EnSightT::VariableTypeT>& vtypes) const;
bool IsVector (const ArrayT<StringT>& labels, int index, StringT& extension, int dof) const;

private:
bool fBinary;
int  fNumDigits;
AutoArrayT<double> fTimeValues;
};

#endif

