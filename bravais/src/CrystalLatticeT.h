/* $Id: CrystalLatticeT.h,v 1.6 2002-08-02 02:07:49 saubry Exp $ */

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
	dArray2DT matrix_rotation;
	double angle_rotation;
	dArrayT norm_vec;
	double density;

public:

	// Constructor 
	CrystalLatticeT(int nlsd, int nuca,
			dArray2DT mat_rot,double angle);
	// Copy Constructor 
	CrystalLatticeT(const CrystalLatticeT& source);
	// Destructor
	~CrystalLatticeT() { }

	int GetNLSD() { return nLSD; }
	int GetNUCA() { return nUCA; }

	double GetAngleRotation() { return angle_rotation; };
	dArray2DT GetMatrixRotation() { return matrix_rotation;};

	virtual const dArrayT& GetLatticeParameters() = 0;
	virtual const dArray2DT& GetBasis() = 0;
	

	void CalculateDensity();
	double GetDensity();


	dArray2DT BasisRotation(dArray2DT A);
	//dArrayT AxisRotation(dArrayT vectorIN);

};

#endif
