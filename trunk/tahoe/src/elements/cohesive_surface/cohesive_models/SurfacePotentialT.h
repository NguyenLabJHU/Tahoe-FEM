/* $Id: SurfacePotentialT.h,v 1.30 2007-03-12 02:20:07 cjkimme Exp $ */
/* created: paklein (06/20/1999) */
#ifndef _SURFACE_POTENTIAL_T_H_
#define _SURFACE_POTENTIAL_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dMatrixT.h"

#include "GlobalT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class StringT;

/** base class for surface potential. The traction evolves in the local
 * frame as a function of the current opening displacement, the rate of
 * opening displacement, and a vector state variables. The "needs" of the
 * model are probed by the host code. The state variables are integrated
 * over the current time step during the call to SurfacePotentialT::Traction. 
 * during all other calls, the state variable array is not modified. */
class SurfacePotentialT: public ParameterInterfaceT
{
public:

	/** surface potential types - derived classes */
	enum CodeT {kXuNeedleman = 0, /**< elastic potential developed by Xu and Needleman */
	    kTvergaardHutchinson = 1, /**< tri-linear potential */
	           kLinearDamage = 2, /**< irreversible linear decay */
	kViscTvergaardHutchinson = 3, /**< T-H with viscous dissipation */
	               kTijssens = 4, /**< Tjissens rate dependent model */
	                kRateDep = 5, /**< simplified rate model */
              kTiedPotential = 6, /**< wrapper for models w/o initial load-up */
              	  kYoonAllen = 7, /**< Allen history-dependent law */
			 kSimoViscoElast = 8, /**< Simo's thermo-visco-elasto model */
           kInelasticDuctile = 9, /**< rate-based ductile fracture model */
           				kMR = 10, /**< Elastoplastic Cohesive Model for Geomaterials*/
           			 kMR_RP = 11, /**< Rigid-plastic Cohesive Model for Geomaterials*/
           		kMR_NodalRP = 12, /**< Nodal-Rigid-plastic Cohesive Model for Geomaterials*/
  kTvergaardHutchinsonRigid = 101, /**< tri-linear potential */
       kInelasticDuctile_RP = 109, /**< rate-based ductile fracture model */
 kTvergaardHutchinsonIrreversible = 110 /**< tri-linear potential with unloading straight to origin to break reversability */
  };

	/** surface element status codes */
	enum StatusT {Precritical = 0, /**< loading phase */
	                 Critical = 1, /**< unloading phase */
	                   Failed = 2  /**< beyond zero-traction opening */};

	/** factory method. Construct a new instance of a sub-class of SurfacePotentialT
	 * with the given ParameterInterfaceT name, or return NULL if the name is
	 * not recognized. */
	static SurfacePotentialT* New(const char* name);

	/** constructor */
	SurfacePotentialT(int ndof);

	/** destructor */
	virtual ~SurfacePotentialT(void);

	/** set the source of the time step */
	virtual void SetTimeStep(const double&) {};

	/** return the number of state variables needed by the model */
	virtual int NumStateVariables(void) const = 0;

	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void InitStateVariables(ArrayT<double>& state);
	
	/** initialize the state variable array. By default, initialization
	 * involves only setting the array to zero. */
	virtual void UpdateStateVariables(ArrayT<double>& state) {/*do nothing by default*/};

	/** dissipated energy. Total amount of energy dissipated reaching
	 * the current state. */
	virtual double FractureEnergy(const ArrayT<double>& state) = 0;

	/*cohesive strength*/
	virtual const double FractureStrength(void) const;

	/** incremental heat. The amount of energy per unit undeformed area
	 * released as heat over the current increment */
	virtual double IncrementalHeat(const dArrayT& jump, const ArrayT<double>& state);

	/** potential energy */
	virtual double Potential(const dArrayT& jump, const ArrayT<double>& state) = 0;
	
	/** surface traction. Internal variables are integrated over the current
	 * time step. */	
	virtual const dArrayT& Traction(const dArrayT& jump, ArrayT<double>& state, const dArrayT& sigma, bool qIntegrate) = 0;

	/** tangent stiffness */
	virtual const dMatrixT& Stiffness(const dArrayT& jump, const ArrayT<double>& state, const dArrayT& sigma) = 0;

	/** type of stiffness matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const { return GlobalT::kSymmetric; }

	/** surface status */
	virtual StatusT Status(const dArrayT& jump, const ArrayT<double>& state) = 0;

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
	//virtual bool NeedsNodalInfo(void);
	//virtual int NodalQuantityNeeded(void);

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

inline const double SurfacePotentialT::FractureStrength(void) const
{
	cout <<"\nSurfacePotentialT::FractureStrength:  Cohesive energy not yet defined.";
	return(0.0);
}

} // namespace Tahoe 
#endif /* _SURFACE_POTENTIAL_T_H_ */
