/* $Id: IsoVIB2D.h,v 1.9 2004-07-15 08:27:45 paklein Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ISO_VIB_2D_H_
#define _ISO_VIB_2D_H_

/* base classes */
#include "FSSolidMatT.h"
#include "VIB.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class IsoVIB2D: public FSSolidMatT, public VIB
{
public:

	/* constructor */
	IsoVIB2D(void);

	/* destructor */
	~IsoVIB2D(void);
	
	/** \name spatial description */
	/*@{*/
	/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void);

	/** Cauchy stress */
	virtual const dSymMatrixT& s_ij(void);

	/** return the pressure associated with the last call to 
	 * SolidMaterialT::s_ij. See SolidMaterialT::Pressure
	 * for more information. */
	virtual double Pressure(void) const { return fEigs.Sum()/3.0; };
	/*@}*/

	/* material description */
	virtual const dMatrixT& C_IJKL(void); // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
//TEMP - not yet optimized for total Lagrangian formulation.
//       calls to these write error message and throw ExceptionT::xception

	/* strain energy density */
	virtual double StrainEnergyDensity(void);

	//TEMP
	const dSymMatrixT& CurvatureTensor(void);

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

protected:

	/* strained lengths in terms of the Lagrangian stretch eigenvalues */
	void ComputeLengths(const dArrayT& eigs);

private:

	/* initialize angle tables */
	void Construct(void);

protected:
	
	/* integration point generator */
	CirclePointsT*	fCircle;
	  	
	/* eigenvalues */
	dArrayT	    fEigs;    // rank2 eigenvalues
	dSymMatrixT	fEigmods; // rank4 eigenvalues
	
	/* spectral decomp solver */
	SpectralDecompT	fSpectral;
	
private:

	/* stretch */
	dSymMatrixT fb;

	/* return values */
	dMatrixT    fModulus;
	dSymMatrixT fStress;
};

} // namespace Tahoe 
#endif /* _ISO_VIB_2D_H_ */
