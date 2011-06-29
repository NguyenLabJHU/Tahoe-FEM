/* $Id: OgdenIsotropicT.h,v 1.9 2004-07-15 08:27:22 paklein Exp $ */
/* created: paklein (10/01/2000) */
#ifndef _OGDEN_ISOTROPIC_T_H_
#define _OGDEN_ISOTROPIC_T_H_

/* base classes */
#include "FSIsotropicMatT.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

/* forward declarations */
class SpectralDecompT;

/** base class for large deformation isotropic material following
 * Ogden's spectral formulation. Derived types need only to overload
 * OgdenIsotropicT::dWdE and OgdenIsotropicT::dWdE. */
class OgdenIsotropicT: public FSIsotropicMatT
{
public:

	/** constructor */
	OgdenIsotropicT(void);

	/** destructor */
	virtual ~OgdenIsotropicT(void);

	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const;
	/*@}*/

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/* principal values of the PK2 stress given principal values of the stretch 
	 * tensors, i.e., the principal stretches squared */
	virtual void dWdE(const dArrayT& eigenstretch2, dArrayT& eigenstress) = 0;
	virtual void ddWddE(const dArrayT& eigenstretch2, dArrayT& eigenstress,
		dSymMatrixT& eigenmod) = 0;

	/* return true of model is purely 2D, plain stress */
	virtual bool PurePlaneStress(void) const { return false; };

private:

	/* construct symmetric rank-4 mixed-direction tensor (6.1.44) */
	void MixedRank4_2D(const dArrayT& a, const dArrayT& b,
		dMatrixT& rank4_ab) const;
	void MixedRank4_3D(const dArrayT& a, const dArrayT& b,
		dMatrixT& rank4_ab) const;

protected:

	/* spectral operations */
	SpectralDecompT* fSpectralDecomp;

	/* work space */
	dSymMatrixT fC;
	dArrayT     fEigs; //TEMP - need this??
	dArrayT     fdWdE;
	dSymMatrixT fddWddE;
	dMatrixT    fModMat;
	
	/* return values */
	dMatrixT    fModulus;
	dSymMatrixT fStress;
};

} // namespace Tahoe 
#endif /* _OGDEN_ISOTROPIC_T_H_ */
