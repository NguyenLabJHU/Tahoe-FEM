/* $Id: TiedPotentialT.h,v 1.17 2003-06-09 06:44:36 paklein Exp $ */
/* created: cjkimme (04/15/2002) */

#ifndef _TIED_POTENTIAL_T_H_
#define _TIED_POTENTIAL_T_H_

/* base classes */
#include "SurfacePotentialT.h"
#include "TiedPotentialBaseT.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class dArray2DT;

/** cohesive potential from Tvergaard and Hutchinson. This model is
 * described in JMPS v41, n6, 1995, 1119-1135. See SurfacePotentialT
 * for more information about the */
class TiedPotentialT: public SurfacePotentialT, public TiedPotentialBaseT
{
public:
	
	enum sbntmaT {kAverageCode = 2};

	/** constructor */
	TiedPotentialT(ifstreamT& in);

	/** return the number of state variables needed by the model */
	int NumStateVariables(void) const;
	
	/** initialize the state variable array. */
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
	 * allocated by the function. Returns empty by default. */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute the output variables.
	 * \param destination of output values. Allocated by the host code */
	virtual void ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
		dArrayT& output);

	virtual bool NeedsNodalInfo(void) const;
	
	virtual int NodalQuantityNeeded(void) const;

	/** rotate nodal values to local frame */
	virtual bool RotateNodalQuantity(void) const { return true; };
	
	virtual bool InitiationQ(const nArrayT<double>& sigma) const;
	
	/** whether or not potential may retie nodes */
	virtual bool NodesMayRetie(void) const;
	
	/** returns true if criterium for retieing is met */
	virtual bool RetieQ(const nArrayT<double>& sigma, const ArrayT<double>& state, const dArrayT& jump_u) const;

	/** location in state variable array of the state flag */
	virtual int TiedStatusPosition(void) const;
	
protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output - FALSE by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;
		
private:

	/* traction potential parameters */
	int qTv;
	
	double d_n;   // characteristic normal opening
	double d_t;   // characteristic tangent opening  	
	double phi_n; // mode I work to fracture
	double r_fail; 
	double fsigma, fL_0, fL_1, fL_2;
	double q, r;
	
	double fnvec1, fnvec2; /*components of direction
	  in which to sample the stress for freeing nodes */
	double fsigma_critical; /* Initiation traction */
	
	bool qRetieNodes; // true if nodes might be retied
};

} // namespace Tahoe 
#endif /* _TIED_POTENTIAL_T_H_ */
