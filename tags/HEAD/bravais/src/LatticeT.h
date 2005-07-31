/* $Id: LatticeT.h,v 1.1.1.1 2002-02-28 02:13:08 jzimmer Exp $ */
#include "StringT.h"
#include "ExceptionCodes.h"
#include "dArrayT.h"
#include "dArray2DT.h"

class ifstreamT;

class LatticeT {
protected:
	int nLSD, nUCA;
	dArray2DT vBasis;
	dArrayT vLatticeParameters;
	dArray2DT mRotation;	// should this be lattice object, 
				// or part of the Box.FillVolume() function?
	double density;
public:
	LatticeT(int nlsd, int nuca);
	~LatticeT() { }

	void Dimension();
	virtual void SetBasis() = 0;
	virtual void SetLatticeParameters(ifstreamT& in) = 0;
	void CalculateDensity();
	double GetDensity();
};

class FCCT : public LatticeT {
public:
	FCCT();
	~FCCT() { }

	void SetBasis();
	void SetLatticeParameters(ifstreamT& in);
};
