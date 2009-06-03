/* $Id: EAMFCC3DMatT_edge.h,v 1.2 2009-06-03 03:20:31 hspark Exp $ */
/* created: hspark (6/2/2009) */
#ifndef _EAMFCC3DMatT_EDGE_H_
#define _EAMFCC3DMatT_EDGE_H_

/* base class */
#include "NL_E_MatT.h"

namespace Tahoe {

/* forward declarations */
class EAMFCC3D_edge;

/** plane strain EAM material */
class EAMFCC3DMatT_edge: public NL_E_MatT
{
public:

	/* constructor */
	EAMFCC3DMatT_edge(void);

	/* destructor */
	virtual ~EAMFCC3DMatT_edge(void);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** thickness of surface layer to subtract off of bulk */
	double SurfaceThickness(void) const { return fSurfaceThickness; };

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

//TEMP
/** spatial tangent modulus */
	virtual const dMatrixT& c_ijkl(void) { return FSSolidMatT::c_ijkl(); };

protected:

	/* compute the symetric Cij reduced index matrix */
	virtual void ComputeModuli(const dSymMatrixT& E, dMatrixT& moduli);	
	
	/* compute the symetric 2nd Piola-Kirchhoff reduced index vector */
	virtual void ComputePK2(const dSymMatrixT& E, dSymMatrixT& PK2);

	/* returns the strain energy density for the specified strain */
	virtual double ComputeEnergyDensity(const dSymMatrixT& E);
	
protected:

	/** surface thickness accesssor */
	double fSurfaceThickness;
	
	/** strain-independent subtraction parameter */
	double fAlpha;
	
	/** Cauchy-Born EAM solver */
	EAMFCC3D_edge* fEAM;
};

} /* namespace Tahoe */

#endif /* _EAMFCC3DMatT_EDGE_H_ */
