/* $Id: VolumeT.h,v 1.2 2002-03-06 01:55:43 jzimmer Exp $ */

#ifndef _VOLUME_T_H_
#define _VOLUME_T_H_

#include <iostream>
#include "dArrayT.h"
#include "dArray2DT.h"
#include "CrystalLatticeT.h"

class ifstreamT;

class VolumeT {
protected:
	int nSD;
	double volume;
	int nATOMS;
	dArray2DT atom_coord;
public:
        VolumeT(int n);
	~VolumeT();
	int GetDimensions();
	double GetVolume();
	void WriteFile();

	bool dimunits;

	virtual void SetSize(ifstreamT& in) = 0;
	virtual void DefineBoundary(CrystalLatticeT* pcl) = 0;
	virtual void CalculateVolume() = 0;
	virtual void FillVolume(CrystalLatticeT* pcl) = 0;
};

#endif

