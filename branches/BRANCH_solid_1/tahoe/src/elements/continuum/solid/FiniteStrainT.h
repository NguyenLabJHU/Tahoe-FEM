/* $Id: FiniteStrainT.h,v 1.1.2.3 2001-06-28 01:24:11 paklein Exp $ */

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

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);

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
  
  	/** return values */
  	ArrayT<dMatrixT> fF_List;
  	ArrayT<dMatrixT> fF_last_List;
};

/* inlines */

inline const dMatrixT& FiniteStrainT::DeformationGradient(void) const
{
	return fF_List[CurrIP()];
}

inline const dMatrixT& FiniteStrainT::DeformationGradient(int ip) const
{
	return fF_List[ip];
}

inline const dMatrixT& FiniteStrainT::DeformationGradient_last(void) const
{
	return fF_last_List[CurrIP()];
}

inline const dMatrixT& FiniteStrainT::DeformationGradient_last(int ip) const
{
	return fF_last_List[ip];
}

#endif /* _FINITE_STRAIN_T_H_ */
