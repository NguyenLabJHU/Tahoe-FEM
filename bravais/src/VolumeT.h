// DEVELOPMENT
/* $Id: VolumeT.h,v 1.9 2003-06-06 23:11:36 saubry Exp $ */

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

	StringT atom_names;
	iArrayT atom_ID;
	dArray2DT atom_coord;
	iArray2DT atom_connectivities;
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

	virtual void CalculateBounds(iArrayT per,CrystalLatticeT* pcl) = 0;
	virtual void CalculateType()=0;
	virtual void CalculatePart()=0;

	virtual iArrayT GetNCells() = 0;
        virtual dArray2DT GetLength() = 0; 

	StringT*   GetAtomNames();
	iArrayT*   GetAtomID();
	dArray2DT* GetAtomCoordinates();
	iArray2DT* GetAtomConnectivities();
	dArray2DT* GetAtomBounds();
	iArrayT*   GetAtomTypes();
};

}
#endif

