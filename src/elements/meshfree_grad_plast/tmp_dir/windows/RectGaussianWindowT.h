/* $Id: RectGaussianWindowT.h,v 1.2 2004-08-25 00:34:39 kyonten Exp $ */

#ifndef _RECT_GAUSSIAN_WINDOW_T_H_
#define _RECT_GAUSSIAN_WINDOW_T_H_

/* base class */
#include "WindowT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"  //kyonten

namespace Tahoe {

/** Rectangular Gaussian window function */
class RectGaussianWindowT: public WindowT
{
  public:

	/** constructor */
	RectGaussianWindowT(const dArrayT& dilation_scaling, double sharpening_factor,
				   double cut_off_factor);

	/** neighbor search type.
	 * \return search type recommended for construction support size */
	virtual SearchTypeT SearchType(void) const { return kConnectivity; };
	
	/** window function name */
	virtual const char* Name(void) const { return "Rectangular Gaussian"; };

	/** shared parameters
	 * \return one parameter for each spatial dimension */
	virtual int NumberOfSupportParameters(void) const { return fDilationScaling.Length(); };

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
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dMatrixT& DDDw); // kyonten (DDDw)

	/* multiple point evaluations */
	virtual int Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw); // kyonten (DDDw)

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
	virtual const dArrayT& RectangularSupportSize(const dArrayT& param_n) const;

	/** spherical support sizes in batch */
	virtual void SphericalSupportSize(const dArray2DT& param_n, ArrayT<double>& support_size) const;

	/** rectangular support sizes in batch */
	virtual void RectangularSupportSize(const dArray2DT& param_n, dArray2DT& support_size) const;
	/*@}*/
	
  private:
  
  	/* window function adjustable parameters */
  	dArrayT fDilationScaling;
  	double fSharpeningFactor;
	double fCutOffFactor;
	
	/* work space */
	dArrayT fNSD;
	dSymMatrixT fNSDsym;
	dMatrixT    fNSDunsym; //for DDDw??
};

} // namespace Tahoe 
#endif /* _RECT_GAUSSIAN_WINDOW_T_H_ */
