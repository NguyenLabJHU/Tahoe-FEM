/* $Id: GaussianWindowT.h,v 1.1 2001-06-13 20:53:57 paklein Exp $ */

#ifndef _GAUSSIAN_WINDOW_T_H_
#define _GAUSSIAN_WINDOW_T_H_

/* base class */
#include "WindowT.h"

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

	/* single point evaluations */
	
	//etc...
	
  private:
  
  	/* window function adjustable parameters */
  	double fDilationScaling;
  	double fSharpeningFudgeFactor;
};

#endif /* _WINDOW_T_H_ */
