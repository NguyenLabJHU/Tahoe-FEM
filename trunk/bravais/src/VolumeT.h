// DEVELOPMENT
/* $Id: VolumeT.h,v 1.16 2003-08-04 23:17:36 saubry Exp $ */

#ifndef _VOLUME_T_H_
#define _VOLUME_T_H_

#include <iostream>

#include "iArrayT.h"
#include "iArray2DT.h"

#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"
#include "ifstreamT.h"

namespace Tahoe {


class VolumeT {

protected:

	int nSD;
	int nATOMS;
	double volume;
	StringT sLATTYPE;

	StringT atom_names;
	iArrayT atom_number;
	ArrayT< StringT > atom_ID;
	dArray2DT atom_coord;
	ArrayT< const iArray2DT * >  atom_connect;
	dArray2DT atom_bounds;
	iArrayT atom_types;
	iArrayT atom_parts;

	iArrayT WhichSort;
	iArrayT Map;

public:

        VolumeT(int n);
	~VolumeT() {};
        VolumeT(const VolumeT& source);

	int GetDimensions();
	int GetNumberAtoms();
	double GetVolume();

	virtual void CreateLattice(CrystalLatticeT* pcl) = 0;
	virtual void SortLattice(CrystalLatticeT* pcl) = 0;

	virtual void CalculateBounds() = 0;

	virtual iArrayT GetNCells() = 0;
        virtual dArray2DT GetLength() = 0; 

	int GetNumberOfAtoms() {return nATOMS;}

	StringT*   GetAtomNames();
	const ArrayT< StringT > *   GetAtomID();
	dArray2DT* GetAtomCoordinates();
	const ArrayT< const iArray2DT * > * GetAtomConnect();
	dArray2DT* GetAtomBounds();
	iArrayT*   GetAtomNumber();
	iArrayT*   GetAtomTypes();
	iArrayT*   GetAtomParts();
};

}
#endif

