/* $Id: XuNeedleman2DT.h,v 1.10.42.2 2004-07-13 16:42:27 paklein Exp $ */
/* created: paklein (11/14/1997) */

#ifndef _XU_NEEDLE_2D_T_H_
#define _XU_NEEDLE_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** Xu-Needleman 2D cohesive surface potential */
class XuNeedleman2DT: public SurfacePotentialT
{
public:

	/** constructor */
	XuNeedleman2DT(void);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const { return 0; };

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate);

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump_u, const ArrayT<double>& state, const dArrayT& sigma);

	/** surface status */
	virtual StatusT Status(const dArrayT& jump_u, const ArrayT<double>& state);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters  */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:

	/** \name traction potential parameters */
	/*@{*/
	double q;     /**< \f$ \phi_t/\phi_n \f$ */
	double r;     /**< \f$ \delta_n^* /\delta_n \f$ */
	
	double d_n;   /**< characteristic normal opening \f$ \delta_n \f$ */
	double d_t;   /**< characteristic tangent opening \f$ \delta_t \f$ */	
	double phi_n; /**< mode I work to fracture \f$ \phi_n \f$ */

	double r_fail; /**< \f$ \Delta/\delta_{n,t} \f$ for which surface is considered failed */
	/*@}*/

	/** \name additional penetration stiffness */
	/*@{*/
	double fKratio; /**< stiffening ratio */
	double fK;
	/*@}*/
};

} // namespace Tahoe 
#endif /* _XU_NEEDLE_2D_T_H_ */
