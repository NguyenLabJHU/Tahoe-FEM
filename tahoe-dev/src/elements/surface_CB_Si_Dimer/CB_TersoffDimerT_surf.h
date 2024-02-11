/* $Id: CB_TersoffDimerT_surf.h,v 1.1 2007-11-09 16:53:55 hspark Exp $ */
/* created: paklein (10/14/1998) */
#ifndef _CB_TERSOFFDIMER_T_SURF_H_
#define _CB_TERSOFFDIMER_T_SURF_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class TersoffDimerSolverT_surf;

class CB_TersoffDimerT_surf: public NL_E_MatT
{
public:

	/** constructor */
	CB_TersoffDimerT_surf(void);

	/* destructor */
	virtual ~CB_TersoffDimerT_surf(void);

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
	
private:

	/* compute the 3D stretch tensor from the 2D reduced index
	 * strain vector (assuming plane strain) */
	void StrainToStretch(const dSymMatrixT& E, dMatrixT& C);
	
private:
	
	/** TersoffDimer surface solver */
	TersoffDimerSolverT_surf* fTersoffDimerSolver_surf;
	
	/* work space */
	dArrayT	 fXsi; //internal DOF vector
	dMatrixT fC;	
	dMatrixT fPK2;		
	
	double fSurfaceThickness;
};

} // namespace Tahoe 
#endif /* _CB_TERSOFFDIMER_T_SURF_H_ */
