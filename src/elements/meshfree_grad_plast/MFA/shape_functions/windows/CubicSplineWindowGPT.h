
#ifndef _CUBIC_SPLINE_WINDOW_GP_T_H_
#define _CUBIC_SPLINE_WINDOW_GP_T_H_

/* base class */
#include "WindowGPT.h"

/* direct members */
#include "dArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h" //kyonten
#include "dSymMatrixT.h"

namespace Tahoe 
{

/** Spherical/circular cubic spline window function */
class CubicSplineWindowGPT: public WindowGPT
{
   public:
   
   /* constructor */
	CubicSplineWindowGPT(double dilation_scaling);
	
	/** window function name */
	virtual const char* Name(void) const { return "Cubic Spline"; };

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
  	double fDilationScaling;
  	double fCutOffFactor;
  	
	/* work space */
	dArrayT     fNSD;
	dSymMatrixT fNSDsym; 
	dMatrixT    fNSDunsym; // for DDDw??

};

} // namespace Tahoe 
#endif /* _CUBIC_SPLINE_WINDOW_GP_T_H_ */
