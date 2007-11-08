/* $Id: CB_WurtziteT.h,v 1.1 2007-11-08 19:32:49 hspark Exp $ */
/* created: paklein (10/14/1998) */
#ifndef _CB_WURTZITE_T_H
#define _CB_WURTZITE_T_H

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class WurtziteSolverT;

class CB_WurtziteT: public NL_E_MatT
{
public:

	/** constructor */
	CB_WurtziteT(void);

	/* destructor */
	virtual ~CB_WurtziteT(void);

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
	
	/** Wurtzite solver */
	WurtziteSolverT* fWurtziteSolver;
	
	/* work space */
	dArrayT	 fXsi; //internal DOF vector
	dMatrixT fC;	
	dMatrixT fPK2;		
};

} // namespace Tahoe 
#endif /* _CB_WURTZITE_T_H_ */
