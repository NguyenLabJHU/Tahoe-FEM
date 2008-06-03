/* $Id: CB_TersoffT_surf.h,v 1.3 2008-06-03 21:30:17 hspark Exp $ */
/* created: paklein (10/14/1998) */
#ifndef _CB_TERSOFF_T_SURF_H_
#define _CB_TERSOFF_T_SURF_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class TersoffSolverT_surf;

class CB_TersoffT_surf: public NL_E_MatT
{
public:

	/** constructor */
	CB_TersoffT_surf(void);

	/* destructor */
	virtual ~CB_TersoffT_surf(void);

	/** \name material output variables */
	/*@{*/
	/** return the number of constitutive model output parameters */
	virtual int NumOutputVariables(void) const;

	/** return the labels for model output parameters */
	virtual void OutputLabels(Tahoe::ArrayT<StringT>& labels) const;

	/** return material output variables */
	virtual void ComputeOutput(Tahoe::dArrayT& output);
	/*@}*/	

	/** \name tangent moduli
	 * TEMP: use finite difference approximation for now.
	 */
	/*@{*/
//	virtual const dMatrixT& c_ijkl(void) { return FSSolidMatT::c_ijkl(); };
//	virtual const dMatrixT& C_IJKL(void) { return FSSolidMatT::C_IJKL(); };
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
 	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

	/** thickness of surface layer to subtract off of bulk */
	double SurfaceThickness(void) const { return fSurfaceThickness; };

protected:

	/* compute the symmetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symmetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
	
	/** strain-independent subtraction parameters */
	double fAlpha;
	double fBeta;
	
	/** zero strain surface stiffness */
	dMatrixT fSS0;
	
private:

	/* compute the 3D stretch tensor from the 2D reduced index
	 * strain vector (assuming plane strain) */
	void StrainToStretch(const dSymMatrixT& E, dMatrixT& C);
	
private:
	
	/** Tersoff surface solver */
	TersoffSolverT_surf* fTersoffSolver_surf;
	
	/* work space */
	dArrayT	 fXsi; //internal DOF vector
	dMatrixT fC;	
	dMatrixT fPK2;		
	
	double fSurfaceThickness;
};

} // namespace Tahoe 
#endif /* _CB_TERSOFF_T_SURF_H_ */
