/* $Id: SurfacePotentialT.h,v 1.8 2001-11-16 00:22:42 cjkimme Exp $ */
/* created: paklein (06/20/1999) */

#ifndef _SURFACE_POTENTIAL_T_H_
#define _SURFACE_POTENTIAL_T_H_

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

/* forward declarations */
#include "ios_fwd_decl.h"
class StringT;

/** base class for surface potential. The traction evolves in the local
 * frame as a function of the current opening displacement, the rate of
 * opening displacement, and a vector state variables. The "needs" of the
 * model are probed by the host code. The state variables are integrated
 * over the current time step during the call to SurfacePotentialT::Traction. 
 * during all other calls, the state variable array is not modified. */
class SurfacePotentialT
{
public:

	/** surface potential types - derived classes */
	enum CodeT {kXuNeedleman = 0, /**< elastic potential developed by Xu and Needleman */
	    kTvergaardHutchinson = 1, /**< tri-linear potential */
	           kLinearDamage = 2, /**< irreversible linear decay */
	kViscTvergaardHutchinson = 3, /**< T-H with viscous dissipation */
	               kTijssens = 4  /**< */};

	/** surface element status codes */
	enum StatusT {Precritical = 0, /**< loading phase */
	                 Critical = 1, /**< unloading phase */
	                   Failed = 2  /**< beyond zero-traction opening */};

	/** constructor */
	SurfacePotentialT(int ndof);

	/** destructor */
	virtual ~SurfacePotentialT(void);

	/** return the number of state variables needed by the model */
	virtual int NumStateVariables(void) const = 0;

	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void InitStateVariables(ArrayT<double>& state);
	
	/** dissipated energy */
	virtual double FractureEnergy(const ArrayT<double>& state) = 0;

	/** potential energy */
	virtual double Potential(const dArrayT& jump, const ArrayT<double>& state) = 0;
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump, ArrayT<double>& state) = 0;

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump, const ArrayT<double>& state) = 0;

	/** surface status */
	virtual StatusT Status(const dArrayT& jump, const ArrayT<double>& state) = 0;
	
	/** write model name to output */
	virtual void PrintName(ostream& out) const = 0;

	/** write model parameters */
	virtual void Print(ostream& out) const = 0;

	/** return the number of output variables */
	virtual int NumOutputVariables(void) const;

	/** return labels for the output variables.
	 * \param labels returns with the labels for the output variables. Space is
	 *        allocate by the function. Returns empty by default. */
	virtual void OutputLabels(ArrayT<StringT>& labels) const;

	/** compute the output variables.
	 * \param destination of output values. Allocated by the host code */
	virtual void ComputeOutput(const dArrayT& jump, const ArrayT<double>& state, 
		dArrayT& output);

	/** returns true if the potential needs access to physical quantities
at the nodes. Returns false by default. */
	virtual bool NeedsNodalInfo(void);
	virtual int NodalQuantityNeeded(void);
	virtual double ComputeNodalValue(const dArrayT &); 
        virtual void UpdateStateVariables(const dArrayT & IPdata, ArrayT<double> &);
	virtual int ElementGroupNeeded(void);

	/** returns true if two materials have compatible nodal outputs */
	static bool CompatibleOutput(const SurfacePotentialT&, const SurfacePotentialT&);

protected:

	/** return true if the potential has compatible (type and sequence)
	 * nodal output, returns false by default */
	virtual bool CompatibleOutput(const SurfacePotentialT& potential) const;

protected:

	/* return values */
	dArrayT  fTraction;  /**< traction return value */
	dMatrixT fStiffness; /**< stiffness return value */
};

#endif /* _SURFACE_POTENTIAL_T_H_ */
