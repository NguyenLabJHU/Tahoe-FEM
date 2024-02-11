/* $Id: FSFiberMatSplitT.h,v 1.1 2010-06-24 14:35:49 thao Exp $ */
/* created: paklein (06/09/1997) */
#ifndef _FD_FIB_MATSplit_T_H_
#define _FD_FIB_MATSplit_T_H_

/* base class */
#include "FSFiberMatT.h"

namespace Tahoe {


/** base class for finite deformation fiber composite constitutive models with deviatoric volumetric split formulation. The interface *
 * provides access to the element-computed fiber orientation vectors in the global (lab) *
 * cartesian coordinates.                                                                */
class FSFiberMatSplitT: public FSFiberMatT
{
public:

	/** constructor */
	FSFiberMatSplitT(void);


	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/** describe the parameters needed by the interface */
	virtual void TakeParameterList(const ParameterListT& list);

protected:

	/*computes eq isotropic matrix stress*/
	virtual void ComputeMatrixStress (const dSymMatrixT& Stretch, dSymMatrixT& Stress);

	/*computes eq matrix moduli*/
	virtual void ComputeMatrixMod (const dSymMatrixT& Stretch, dSymMatrixT& Stress, dMatrixT& Mod);

	/*subsequent derived classes must define the following functions*/
	/*calculates  matrix contribution to 2PK stress*/
	virtual void ComputeDevMatrixStress(const dSymMatrixT& Cbar,  dSymMatrixT& Stress)=0;

	/*calculates matrix contribution to modulus*/
	virtual void ComputeDevMatrixMod(const dSymMatrixT& Cbar, dSymMatrixT& Stress, dMatrixT& Mod)=0;
	
	/*calculates  matrix contribution to 2PK stress*/
	virtual double ComputeVolMatrixStress(const double I3)=0;

	/*calculates matrix contribution to modulus*/
	virtual double ComputeVolMatrixMod(const double I3)=0;


private:
	/* stretch */
	dSymMatrixT fCbar;
	dSymMatrixT fInverse;
	double fI3;
	
	/* return values */
	dMatrixT	fModbar;
	
	dSymMatrixT fSbar;

	/*workspace*/
	dSymMatrixT fSymMat1;
	dSymMatrixT fSymMat2;
	dMatrixT fMat;
	
};

} /* namespace Tahoe */

#endif /* _FD_STRUCT_MAT_T_H_ */
