/* $Id: ModCBSolverT.h,v 1.4 2004-07-15 08:28:36 paklein Exp $ */
/* created: paklein (05/27/1997) */
#ifndef _MODCB_SOLVER_T_H_
#define _MODCB_SOLVER_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "iArray2DT.h"
#include "LengthsAndAnglesT.h"
#include "SWDataT.h" //TEMP

namespace Tahoe {

/* forward declaration */
class ThermalDilatationT;
class TwoBodyT;
class ThreeBodyT;

class ModCBSolverT: public ParameterInterfaceT
{
public:

	/* potential type */
	enum PotentialTypeT {kSW = 0,
                       kPTHT = 1,
                    kTersoff = 2};

	/** constructor */
	ModCBSolverT(const ThermalDilatationT* thermal);

	/** destructor */
	~ModCBSolverT(void);

	/** \name constitutive properties */
	/*@{*/
	/** moduli - C_IJKL */
	void SetModuli(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& moduli);

	/** stress - S_IJ (2nd PK) */
	void SetStress(const dMatrixT& CIJ, dArrayT& Xsi, dMatrixT& stress);

	/** strain energy density */
	double StrainEnergyDensity(const dMatrixT& CIJ, dArrayT& Xsi);
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
 
  	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

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

	/** lattice geometry */
	LengthsAndAnglesT* fGeometry;
	
	/** thermal dilatation */
	const ThermalDilatationT* fThermal;
	
	/* potential functions and derivatives */
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

} /* namespace Tahoe */

#endif /* _MODCB_SOLVER_T_H_ */
