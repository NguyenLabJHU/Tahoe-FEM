// DEVELOPMENT
/* $Id: CrystalLatticeT.h,v 1.10 2003-06-06 23:11:36 saubry Exp $ */

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

	dArray2DT vBasis;            // atoms in cell
	dArrayT vLatticeParameters;  // lattice parameters
	dArray2DT vAxis;             // Bravais vectors

	int WhichRot;
	dArray2DT matrix_rotation;
	double angle_rotation;
	dArrayT norm_vec;

	double density;

public:

	// Constructor 
	CrystalLatticeT(int nlsd, int nuca,int which_rot,
			dArray2DT mat_rot,double angle);
	// Copy Constructor 
	CrystalLatticeT(const CrystalLatticeT& source);
	// Destructor
	~CrystalLatticeT() { }

	int GetNLSD() { return nLSD; }
	int GetNUCA() { return nUCA; }

	int GetRotMeth() { return WhichRot; };
	double GetAngleRotation() { return angle_rotation; };
	dArray2DT GetMatrixRotation() { return matrix_rotation;};

	virtual const dArrayT& GetLatticeParameters() = 0;
	virtual const dArray2DT& GetBasis() = 0;
	virtual const dArray2DT& GetAxis() = 0;
	

	void   CalculateDensity();
	double GetDensity();

	dArray2DT AxisRotation(dArray2DT A);
	dArrayT VectorRotation(dArrayT v);
};

#endif
