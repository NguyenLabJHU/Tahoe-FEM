/* $Id: BoxT.h,v 1.1 2002-03-06 01:55:43 jzimmer Exp $ */

#ifndef _BOX_T_H_
#define _BOX_T_H_

#include <iostream>
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "VolumeT.h"
#include "CrystalLatticeT.h"

class ifstreamT;

class BoxT : public VolumeT {
protected:
	iArrayT ncells;
        dArrayT length;
        dArray2DT surfaces;
public:
        BoxT(int n);
        ~BoxT();
	void SetSize(ifstreamT& in);
        void DefineBoundary(CrystalLatticeT* pcl);
        void CalculateVolume();
        void FillVolume(CrystalLatticeT* pcl);
};



#endif

