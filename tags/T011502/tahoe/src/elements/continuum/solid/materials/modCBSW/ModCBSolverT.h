/* $Id: ModCBSolverT.h,v 1.1.1.1 2001-01-29 08:20:26 paklein Exp $ */
/* created: paklein (05/27/1997)                                          */
/* Q defines the orientation of the crystals' natural coordinates         */
/* and the global coordinate frame. Q is defined as:                      */
/* 			Q = d x_natural / d x_global                                        */
/* So that the vectors are transformed by:                                */
/* 			r_global = Transpose[Q].r_natural                                   */
/* NOTE: All calculations done in 3D.                                     */

#ifndef _MODCB_SOLVER_T_H_
#define _MODCB_SOLVER_T_H_

/* direct members */
#include "iArray2DT.h"
#include "LengthsAndAnglesT.h"
#include "SWDataT.h" //TEMP

/* forward declaration */
class ThermalDilatationT;
class TwoBodyT;
class ThreeBodyT;

class ModCBSolverT
{
public:

	/* potential type */
	enum PotentialTypeT {kSW = 0,
                       kPTHT = 1,
                    kTersoff = 2};

	/* Constructor */
	ModCBSolverT(const dMatrixT& Q, const ThermalDilatationT* thermal,
		ifstreamT& in, bool equilibrate);

	/* Destructor */
	~ModCBSolverT(void);
	
	/* moduli - C_IJKL */
	void SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli);

	/* stress - S_IJ (2nd PK) */
	void SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress);

	/* strain energy density */
	double StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi);

	/*
	 * Printing parameters.
	 */
	void Print(ostream& out) const;
	void PrintName(ostream& out) const;

private:

	/* Minimize the energy wrt Xsi using the initial value passed */
	void Equilibrate(const dMatrixT& CIJ, dArrayT& Xsi);

	/* set free dof - triggers recomputation */
	void SetdXsi(const dMatrixT& CIJ, const dArrayT& Xsi);
	void SetAll(const dMatrixT& CIJ);		

private:

	/* mod CB flag */
	bool fEquilibrate;

	/* pairs for 3-body potentials */
	const iArray2DT fPairs;

	/* lattice geometry */
	LengthsAndAnglesT	fGeometry;
	
	/* potential functions and derivatives */
	int			fPotential;
	SWDataT		fSW;	//should really make class to manage
	TwoBodyT*	f2Body; //2 and 3 body potentials together
	ThreeBodyT*	f3Body;

	/* derivatives wrt. Xsi */
	dArrayT		dXsi;
	dMatrixT	dXsidXsi;
	
	/* derivatives wrt. C */
	dMatrixT	dCdC_hat;

	/* mixed derivatives wrt. C and Xsi */
	dMatrixT	dCdXsi_hat;
	
	/* work space */
	ArrayT<nArrayT<double>*>  fMatrices; //linear combo
	dMatrixT	fMat1, fMat2;
	dMatrixT	fGradl_i;
	dArrayT		fVec;
	dSymMatrixT	fSymMat1;
	dMatrixT	fTempRank4;
	dMatrixT	fTempMixed;
	dMatrixT	fGradl_C;
};

#endif /* _MODCB_SOLVER_T_H_ */
