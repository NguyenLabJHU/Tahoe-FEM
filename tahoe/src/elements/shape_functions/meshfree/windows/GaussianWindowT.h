/* $Id: GaussianWindowT.h,v 1.13 2004-10-31 20:48:41 paklein Exp $ */

#ifndef _GAUSSIAN_WINDOW_T_H_
#define _GAUSSIAN_WINDOW_T_H_

/* base class */
#include "WindowT.h"

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/** Spherical/circular Gaussian window function */
class GaussianWindowT: public WindowT
{
  public:

	/** constructor */
	GaussianWindowT(double dilation_scaling, double sharpening_factor,
		double cut_off_factor);
	
	/** window function name */
	virtual const char* Name(void) const { return "Gaussian"; };

	/** neighbor search type.
	 * \return search type recommended for construction support size */
	virtual SearchTypeT SearchType(void) const { return kSpherical; };

	/** shared parameters
	 * \return 1 since the function varies from node to node only
	 * depending on the support radius. */
	virtual int NumberOfSupportParameters(void) const { return 1; };

	/** "synchronization" of nodal field parameters. 
	 * Take "max" over both sets for each node
	 * \params params_1 support size set 1
	 * \params params_2 support size set 2 */
	virtual void SynchronizeSupportParameters(dArray2DT& params_1, 
		dArray2DT& params_2) const;

	/** write parameters to output stream */
	virtual void WriteParameters(ostream& out) const;

	/* single point evaluations */
	virtual bool Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw);

	/* multiple point evaluations */
	virtual int Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw);

	/** \name coverage tests */
	/*@{*/
	/** single point */
	virtual bool Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const;

	/** multiple points */
	virtual int Covers(const dArray2DT& x_n, const dArrayT& x, 
		const dArray2DT& param_n, ArrayT<bool>& covers) const;
	/*@}*/

	/** support dimensions */
	/*@{*/
	/** spherical upport size */
	virtual double SphericalSupportSize(const dArrayT& param_n) const;

	/** rectangular support size */
	virtual void RectangularSupportSize(const dArrayT& param_n, dArrayT& support_size) const;

	/** spherical support sizes in batch */
	virtual void SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const;
	/*@}*/

  private:
  
  	/** \name window function adjustable parameters */
  	/*@{*/
  	double fDilationScaling;
  	double fSharpeningFactor;
  	double fCutOffFactor;
  	/*@}*/
  	
	/* work space */
	dArrayT     fNSD;
	dSymMatrixT fNSDsym;
	dArrayT fSupportSize;
};

} /* namespace Tahoe */

#endif /* _GAUSSIAN_WINDOW_T_H_ */
