/* $Id: FiniteStrainT.h,v 1.1.2.1 2001-06-22 01:21:21 paklein Exp $ */

#ifndef _FINITE_STRAIN_T_H_
#define _FINITE_STRAIN_T_H_

/* base class */
#include "ElasticT.h"

/** Interface for linear strain deformation and field gradients */
class FiniteStrainT: public ElasticT
{
  public:
      
	/** constructor */
	FiniteStrainT(FEManagerT& fe_manager);

	/** total deformation gradient */
	const dMatrixT& DeformationGradient(void) const;
	const dMatrixT& DeformationGradient(int ip) const;

	/** total strain from the end of the previous time step */
	const dMatrixT& DeformationGradient_last(void) const;
	const dMatrixT& DeformationGradient_last(int ip) const;

	/** compute field gradients with respect to current coordinates */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients with respect to reference coordinates */
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient_reference(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;
};

#endif /* _FINITE_STRAIN_T_H_ */
