/* $Id: FCCT.cpp,v 1.1 2002-03-06 01:55:43 jzimmer Exp $ */
#include "FCCT.h"
#include "CrystalLatticeT.h"
#include <iostream>
#include <fstream>
#include "StringT.h"
#include "ifstreamT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

FCCT::FCCT() : CrystalLatticeT(3,4) {
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
	for (int i=0; i<nLSD; i++) {
		vLatticeParameters[i] = alat;
	}
}

const dArray2DT& FCCT::GetBasis() {
	if (!&vBasis) {
		cout << "vBasis doesn't exist!" << "\n";
		throw eDatabaseFail;
	}
	else {
		return vBasis;		
	}

}
const dArrayT& FCCT::GetLatticeParameters() {
	if (!&vLatticeParameters) {
		cout << "vLatticeParameters doesn't exist!" << "\n";	
		throw eDatabaseFail;
	}
	else {
		return vLatticeParameters;
	}
}	


