/* $Id: TvergHutchRigid2DT.h,v 1.1 2003-08-08 00:43:00 paklein Exp $ */
#ifndef _TVERG_HUTCH_RIGID_2D_T_H_
#define _TVERG_HUTCH_RIGID_2D_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;

/** cohesive potential from Tvergaard and Hutchinson. This model is
 * based described in JMPS v41, n6, 1995, 1119-1135. The model differs
 * in that the parameters have been redfined for an initially rigid response
 * and more flexible shape factor parameters. The model returns the initiation
 * traction, which must be written into the state variables, at zero opening.*/
class TvergHutchRigid2DT: public SurfacePotentialT
{
public:

	/** constructor */
	TvergHutchRigid2DT(ifstreamT& in);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const { return 2; };

	/** initialize the state variable array */
	virtual void InitStateVariables(ArrayT<double>& state);

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

	/** write model name to output */
	virtual void PrintName(ostream& out) const;

	/** write model parameters */
	virtual void Print(ostream& out) const;

	/** return the number of output variables. returns 0 by default. */
	virtual int NumOutputVariables(void) const;

	/** return labels for the output variables.
	 * \param labels returns with the labels for the output variables. Space is
	 *        allocate by the function. Returns empty by default. */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute the output variables.
	 * \param destination of output values. Allocated by the host code */
	virtual void ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
		dArrayT& output);

protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;
	
private:

	/** \name traction potential parameters */
	/*@{*/
	double fsigma_max; /**< cohesive stress */
	double fd_c_n;     /**< characteristic normal opening to failure */
	double fd_c_t;     /**< characteristic tangential opening to failure */
	/*@}*/
	
	/** \name non-dimensional opening parameters */
	/*@{*/
	double fL_1;    /**< shape factor from unmodified T-H model that is used
	                 * to compute the stiffness of the compressive response */

	double fL_2;    /**< shape factor opening */
	double fT_2;    /**< traction at TvergHutchRigid2DT::fL_2 */

	double fL_fail; /**< non-dimensional opening to irreversible failure */
	/*@}*/

	/** \name penetration stiffness */
	/*@{*/
	double fpenalty; /**< stiffening multiplier */
	double fK;       /**< penetration stiffness calculated as a function of penalty
	                  * and the initial stiffness of the cohesive potential */
	/*@}*/

};

} /* namespace Tahoe */

#endif /* _TVERG_HUTCH_RIGID_2D_T_H_ */
