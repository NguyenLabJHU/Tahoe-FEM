/* $Id: IsoVIB2D.h,v 1.8 2003-01-29 07:34:52 paklein Exp $ */
/* created: paklein (11/08/1997) */
#ifndef _ISO_VIB_2D_H_
#define _ISO_VIB_2D_H_

/* base classes */
#include "FSSolidMatT.h"
#include "Material2DT.h"
#include "VIB.h"

/* direct members */
#include "SpectralDecompT.h"

namespace Tahoe {

/* forward declarations */
class CirclePointsT;

/** 2D Isotropic VIB solver using spectral decomposition formulation */
class IsoVIB2D: public FSSolidMatT, public Material2DT, public VIB
{
public:

	/* constructor */
	IsoVIB2D(ifstreamT& in, const FSMatSupportT& support);

	/* destructor */
	~IsoVIB2D(void);
	
	/* print parameters */
	virtual void Print(ostream& out) const;
	virtual void PrintName(ostream& out) const;	

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
