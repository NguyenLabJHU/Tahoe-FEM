// DEVELOPMENT
/* $Id: VolumeT.h,v 1.14 2003-08-01 22:54:39 saubry Exp $ */

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
	iArrayT atom_ID;
	dArray2DT atom_coord;
	iArray2DT atom_connectivities;
	dArray2DT atom_bounds;
	iArrayT atom_types;
	iArrayT atom_parts;

	// try to use the other OutputSeT constructor
	ArrayT< const iArray2DT * > atom_array_connect;   
	ArrayT< StringT > atom_array_ID;

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

	virtual void CalculateBounds(CrystalLatticeT* pcl) = 0;

	virtual iArrayT GetNCells() = 0;
        virtual dArray2DT GetLength() = 0; 

	StringT*   GetAtomNames();
	iArrayT*   GetAtomID();
	dArray2DT* GetAtomCoordinates();
	iArray2DT* GetAtomConnectivities();
	dArray2DT* GetAtomBounds();
	iArrayT*   GetAtomTypes();
	iArrayT*   GetAtomParts();


	const ArrayT< const iArray2DT * > * GetAtomArrayConnect();
	const ArrayT< StringT > * GetAtomArrayID();

};

}
#endif

