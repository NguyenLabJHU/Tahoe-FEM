/* $Id: TiedPotentialT.h,v 1.4 2002-07-05 22:27:59 paklein Exp $ */
/* created: cjkimme (04/15/2002) */

#ifndef _TIED_POTENTIAL_T_H_
#define _TIED_POTENTIAL_T_H_

/* base class */
#include "SurfacePotentialT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArray2DT;

/** cohesive potential from Tvergaard and Hutchinson. This model is
 * described in JMPS v41, n6, 1995, 1119-1135. See SurfacePotentialT
 * for more information about the */
class TiedPotentialT: public SurfacePotentialT
{
public:
	
	enum sbntmaT {kAverageCode = 2};

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
	
	static bool InitiationQ(const double* sigma);
	
protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;
		
private:

	/** the time step */
	const double& fTimeStep;

	/* traction potential parameters */
	double q;     // phi_t/phi_n
	double r;     // delta_n* /d_n
	
	double d_n;   // characteristic normal opening
	double d_t;   // characteristic tangent opening  	
	double phi_n; // mode I work to fracture

	double r_fail; // d/d_(n/t) for which surface is considered failed

/* additional penetration stiffness */
	double fKratio; // stiffening ratio
	double fK;
	double fu_t0, fu_n0; /* Offsets of gap vector when nodes are untied */

	static double fsigma_critical;
};

} // namespace Tahoe 
#endif /* _TIED_POTENTIAL_T_H_ */
