/* $Id: GaussianWindowT.h,v 1.3 2001-06-15 14:45:47 hspark Exp $ */

#ifndef _GAUSSIAN_WINDOW_T_H_
#define _GAUSSIAN_WINDOW_T_H_

/* base class */
#include "WindowT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dSymMatrixT.h"

/** Spherical/circular Gaussian window function */
class GaussianWindowT
{
  public:

	/** constructor */
	GaussianWindowT(double dilation_scaling, double sharpening_factor);
	
	/** window function name */
	virtual const char* Name(void) const { return "Gaussian"; };

	/** shared parameters
	 * \return 1 since the function varies from node to node only
	 * depending on the support radius. */
	virtual int NumberOfParameters(void) const { return 1; };
	
	virtual void WriteParameters(ostream& out) const;

	/* single point evaluations */
	virtual void window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dArrayT& DDw);

	/* multiple point evaluations */
	virtual int window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw);

	/* coverage tests */
	/* single point */
	virtual bool Covers(const dArrayT& x_n, const dArrayT& x);

	/* multiple points */
	virtual void Covers(const dArray2DT& x_n, const dArrayT& x, const ArrayT<bool>& covers);
	
	//etc...
	
  private:
  
  	/* window function adjustable parameters */
  	double fDilationScaling;
  	double fSharpeningFudgeFactor;
};

#endif /* _WINDOW_T_H_ */
