/* $Id: CrystalLatticeT.cpp,v 1.1 2002-03-06 01:55:43 jzimmer Exp $ */
#include "CrystalLatticeT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

CrystalLatticeT::CrystalLatticeT(int nlsd, int nuca) {
	nLSD = nlsd;
	nUCA = nuca;
        vBasis.Dimension(nLSD,nUCA);
        vLatticeParameters.Dimension(nLSD);
}

void CrystalLatticeT::CalculateDensity() {
	double ucvolume = 1.0;
	for (int i=0; i<nLSD; i++) {
		ucvolume *= vLatticeParameters[i];
	}
	density = ((double) nUCA)/ucvolume;
}

double CrystalLatticeT::GetDensity() {
	return density;
}


