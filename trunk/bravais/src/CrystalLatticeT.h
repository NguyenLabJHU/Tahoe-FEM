/* $Id: CrystalLatticeT.h,v 1.3 2002-07-24 23:14:56 saubry Exp $ */

#ifndef _CRYSTAL_LATTICE_T_H_
#define _CRYSTAL_LATTICE_T_H_

#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "ifstreamT.h"

using namespace Tahoe;

class CrystalLatticeT {

protected:

	int nLSD, nUCA;
	dArray2DT vBasis;
	dArrayT vLatticeParameters;
	double density;

public:

	// Constructor 
	CrystalLatticeT(int nlsd, int nuca);
	// Copy Constructor 
	CrystalLatticeT(const CrystalLatticeT& source);
	// Destructor
	~CrystalLatticeT() { }

	int GetNLSD() { return nLSD; }
	int GetNUCA() { return nUCA; }

	virtual const dArrayT& GetLatticeParameters() = 0;
	virtual const dArray2DT& GetBasis() = 0;

	void CalculateDensity();
	double GetDensity();
};

#endif
