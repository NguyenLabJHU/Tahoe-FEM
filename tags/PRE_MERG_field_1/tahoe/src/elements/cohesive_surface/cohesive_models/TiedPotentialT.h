/* $Id: TiedPotentialT.h,v 1.1 2002-04-17 15:40:13 cjkimme Exp $ */
/* created: cjkimme (04/15/2002) */

#ifndef _TIED_POTENTIAL_T_H_
#define _TIED_POTENTIAL_T_H_

/* base class */
#include "SurfacePotentialT.h"

/* forward declarations */
class ifstreamT;

/** cohesive potential from Tvergaard and Hutchinson. This model is
 * described in JMPS v41, n6, 1995, 1119-1135. See SurfacePotentialT
 * for more information about the */
class TiedPotentialT: public SurfacePotentialT
{
public:

	/** constructor */
	TiedPotentialT(ifstreamT& in, const double &fTimeStep);

	virtual void InitStateVariables(ArrayT<double>& state);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;

	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump_u, const ArrayT<double>& state);
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump_u, ArrayT<double>& state, const dArrayT& sigma);

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
	 * allocated by the function. Returns empty by default. */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute the output variables.
	 * \param destination of output values. Allocated by the host code */
	virtual void ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
		dArrayT& output);

	virtual bool NeedsNodalInfo(void);
	virtual int NodalQuantityNeeded(void);
	//        virtual double ComputeNodalValue(const dArrayT &);
	//        virtual void UpdateStateVariables(const dArrayT &, ArrayT<double> &);
	virtual int ElementGroupNeeded(void);

protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;
	
private:

	/* traction potential parameters */
	double fsigma_max; /**< cohesive stress */
	double fd_c_n;     /**< characteristic normal opening to failure */
	double fd_c_t;     /**< characteristic tangential opening to failure */
	
	/* non-dimensional opening parameters */
	double fL_1;    /**< non-dimensional opening to initial peak traction */
	double fL_2;    /**< non-dimensional opening to final peak traction */
	double fL_fail; /**< non-dimensional opening to irreversible failure */

	/* penetration stiffness */
	double fpenalty; /**< stiffening multiplier */
	double fK;       /**< penetration stiffness calculated as a function of penalty
	                  * and the initial stiffness of the cohesive potential */
	const double& fTimeStep;
	int fGroup;
	//bool initiationQ;
	double L_2_b, L_2_m;/* fitting constants for rate dependence of L_2 */
	double fslope; /*slope of the 'plateau' of the Tverg-Hutch potential in units of fsigma_max */
	double fu_t0, fu_n0; /* Offsets of gap vector when nodes are untied */
};

#endif /* _TIED_POTENTIAL_T_H_ */