/* $Id: VolumeT.h,v 1.4 2002-07-24 23:14:56 saubry Exp $ */

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

public:

        VolumeT(int n);
	~VolumeT() {};
        VolumeT(const VolumeT& source);

	int GetDimensions();
	int GetNumberAtoms();

	double GetVolume();

	StringT* GetAtomNames();
	iArrayT* GetAtomID();
	dArray2DT* GetAtomCoordinates();
	iArray2DT* GetAtomConnectivities();

	virtual void CalculateVolume() = 0;
	virtual void CreateLattice(CrystalLatticeT* pcl) = 0;

	virtual iArrayT GetNCells() = 0;
        virtual dArrayT GetLength() = 0; 

};

}
#endif

