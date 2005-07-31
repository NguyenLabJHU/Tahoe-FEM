/* $Id: LatticeT.cpp,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $ */
#include "LatticeT.h"

#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

LatticeT::LatticeT(int nlsd, int nuca) {
	nLSD = nlsd;
	nUCA = nuca;
}

void LatticeT::Dimension() {
	vBasis.Dimension(nLSD,nUCA);
	vLatticeParameters.Dimension(nLSD);
}

void LatticeT::CalculateDensity() {
	double ucvolume = 1.0;
	for (int i=0; i<nLSD; i++) {
		ucvolume *= vLatticeParameters[i];
	}
	density = ((double) nUCA)/ucvolume;
}

double LatticeT::GetDensity() {
	return density;
}

FCCT::FCCT() : LatticeT(3,4) {
}

void FCCT::SetBasis() {
	vBasis(0,0) = 0.0;
	vBasis(1,0) = 0.0;
	vBasis(2,0) = 0.0;
	vBasis(0,1) = 0.5;
	vBasis(1,1) = 0.5;
	vBasis(2,1) = 0.0;
	vBasis(0,2) = 0.5;
	vBasis(1,2) = 0.0;
	vBasis(2,2) = 0.5;
	vBasis(0,3) = 0.0;
	vBasis(1,3) = 0.5;
	vBasis(2,3) = 0.5;
}

void FCCT::SetLatticeParameters(ifstreamT& in) {
	double alat;
	in >> alat;
	for (int i=0; i<3; i++) {
		vLatticeParameters[i] = alat;
	}
}
