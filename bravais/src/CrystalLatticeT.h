/* $Id: CrystalLatticeT.h,v 1.1 2002-03-06 01:55:43 jzimmer Exp $ */

#ifndef _CRYSTAL_LATTICE_T_H_
#define _CRYSTAL_LATTICE_T_H_

#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

class ifstreamT;

class CrystalLatticeT {
protected:
	int nLSD, nUCA;
	dArray2DT vBasis;
	dArrayT vLatticeParameters;
	dArray2DT mRotation;	// should this be lattice object, 
				// or part of the Box.FillVolume() function?
	double density;
public:
	CrystalLatticeT(int nlsd, int nuca);
	~CrystalLatticeT() { }

	int GetNLSD() { return nLSD; }
	int GetNUCA() { return nUCA; }
	virtual void SetBasis() = 0;
	virtual void SetLatticeParameters(ifstreamT& in) = 0;
	virtual const dArrayT& GetLatticeParameters() = 0;
	virtual const dArray2DT& GetBasis() = 0;
	void CalculateDensity();
	double GetDensity();
};

#endif
