/* $Id: FiniteStrainT.h,v 1.1.2.2 2001-06-26 07:17:33 paklein Exp $ */

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

  protected:

	/** increment current element */
	virtual bool NextElement(void);	

  private:
  
  	/** matrix indicies */
  	enum MatrixIndexT {kF = 0,
  	                kF_ip = 1,
  	              kF_last = 2,
  	           KF_last_ip = 3};

  private:
  
  	/** return values */
  	ArrayT<dMatrixT> fMatrixList;
  	
  	/** cached matrices flags */
  	iArrayT fIPSet;
};

#endif /* _FINITE_STRAIN_T_H_ */
