/* $Id: SmallStrainT.h,v 1.1.2.2 2001-06-26 07:17:33 paklein Exp $ */

#ifndef _SMALL_STRAIN_T_H_
#define _SMALL_STRAIN_T_H_

/* base class */
#include "ElasticT.h"

/** Interface for linear strain deformation and field gradients */
class SmallStrainT: public ElasticT
{
  public:
      
	/** constructor */
	SmallStrainT(FEManagerT& fe_manager);

	/** total strain */
	const dSymMatrixT& LinearStrain(void) const;
	const dSymMatrixT& LinearStrain(int ip) const;

	/** total strain from the end of the previous time step */
	const dSymMatrixT& LinearStrain_last(void) const;
	const dSymMatrixT& LinearStrain_last(int ip) const;

	/** compute field gradients */
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

	/** compute field gradients from the end of the previous time step */
	void ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u) const;
	void ComputeGradient_last(const LocalArrayT& u, dMatrixT& grad_u, int ip) const;

  protected:

	/** increment current element */
	virtual bool NextElement(void);	

  private:
  
  	/** matrix indicies */
  	enum MatrixIndexT {kstrain = 0,
  	                kstrain_ip = 1,
  	              kstrain_last = 2,
  	           Kstrain_last_ip = 3};
  	           
  private:
  
  	/** return values */
  	ArrayT<dSymMatrixT> fMatrixList;
  	
  	/** cached matrices flags */
  	iArrayT fIPSet;
};

#endif /* _SMALLSTRAIN_T_H_ */
