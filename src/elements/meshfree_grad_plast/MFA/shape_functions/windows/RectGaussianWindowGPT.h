/* $Id: RectGaussianWindowGPT.h,v 1.2 2004-07-14 19:50:31 kyonten Exp $ */

#ifndef _RECT_GAUSSIAN_WINDOW_GP_T_H_
#define _RECT_GAUSSIAN_WINDOW_GP_T_H_

/* base class */
#include "WindowGPT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dSymMatrixT.h"
#include "dMatrixT.h"

namespace Tahoe {

/** Rectangular Gaussian window function */
class RectGaussianWindowGPT: public WindowGPT
{
  public:

	/** constructor */
	RectGaussianWindowGPT(const dArrayT& dilation_scaling, double sharpening_factor,
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

	/** modify nodal shape function parameters */
	virtual void ModifySupportParameters(dArray2DT& nodal_params) const;
	
	/** write parameters to output stream */
	virtual void WriteParameters(ostream& out) const;

	/* single point evaluations */
	virtual bool Window(const dArrayT& x_n, const dArrayT& param_n, const dArrayT& x,
		int order, double& w, dArrayT& Dw, dSymMatrixT& DDw, dMatrixT& DDDw); // kyonten (DDDw)

	/* multiple point evaluations */
	virtual int Window(const dArray2DT& x_n, const dArray2DT& param_n, const dArrayT& x,
		int order, dArrayT& w, dArray2DT& Dw, dArray2DT& DDw, dArray2DT& DDDw); // kyonten (DDDw)

	/* coverage tests */
	
	/* single point */
	virtual bool Covers(const dArrayT& x_n, const dArrayT& x, const dArrayT& param_n) const;

	/* multiple points */
	virtual int Covers(const dArray2DT& x_n, const dArrayT& x, 
		const dArray2DT& param_n, ArrayT<bool>& covers) const;
	
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
#endif /* _RECT_GAUSSIAN_WINDOW_GP_T_H_ */
