/* $Id: FBC_ControllerT.h,v 1.16 2004-09-16 16:49:31 paklein Exp $ */
/* created: paklein (11/17/1997) */
#ifndef _FBC_CONTROLLER_T_H_
#define _FBC_CONTROLLER_T_H_

/* base class */
#include "ParameterInterfaceT.h"

#include "Environment.h"
#include "GlobalT.h"

#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class ifstreamT;
class FEManagerT;
class SolverT;
template <class TYPE> class AutoArrayT;
class iArray2DT;
template <class TYPE> class RaggedArray2DT;
class eIntegratorT;
class StringT;
class FieldT;
class FieldSupportT;

/** base class for all force BC controllers */
class FBC_ControllerT: public ParameterInterfaceT
{
public:

	/** controller codes - derived classes */
	enum CodeT {       kNone =-1,
	            kPenaltyWall = 0,
	          kPenaltySphere = 1,
               kAugLagSphere = 2,
            kMFPenaltySphere = 3,
                 kAugLagWall = 4,
            kPenaltyCylinder = 5,
               kMFAugLagMult = 6,
             kAugLagCylinder = 7};

	/** converts strings to FBC_ControllerT::CodeT */
	static CodeT Code(const char* name);

	/** constructor */
	FBC_ControllerT(void);

	/* destructor */
	virtual ~FBC_ControllerT(void);

	/** set the associated field */
	virtual void SetField(const FieldT& field);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const = 0;

	/* nodally generated DOF's and tags (pseudo nodes) */
	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
	virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
		AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
		AutoArrayT<const iArray2DT*>& equivalent_nodes) const;

	/* initial condition/restart functions
	 *
	 * Set to initial conditions.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	virtual void InitialCondition(void) = 0;
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/** \name apply force and tangent contributions */
	/*@{*/
	/** tangent
	 * \param sys_type "maximum" tangent type needed by the solver. The GlobalT::SystemTypeT
	 *        enum is ordered by generality. The solver should indicate the most general
	 *        system type that is actually needed. */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type) = 0;
	virtual void ApplyRHS(void) = 0;

	/** returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);
	/*@}*/

	/* initialize/finalize step */
	virtual void InitStep(void) = 0;
	virtual void CloseStep(void) = 0;

	/* reset to the last known solution */
	virtual void Reset(void) = 0;

	/** \name writing results */
	/*@{*/
	/** register data for output */
	virtual void RegisterOutput(void) = 0;

	/** write results */
	virtual void WriteOutput(ostream& out) const = 0;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** accept parameter list. Must be called after FBC_ControllerT::SetField */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** return the FieldT or throw exception if not set */
	const FieldT& Field(void) const;

	/** return the FieldSupportT or throw exception if not set */
	const FieldSupportT& FieldSupport(void) const;

protected:

	/** support */
	const FieldSupportT* fFieldSupport;

	/** the field */
	const FieldT* fField;

	/** equation group */
	int fGroup;

	/** element time integration parameters */
	const eIntegratorT* fIntegrator;
};

/* return the FieldT or throw exception if not set */
inline const FieldT& FBC_ControllerT::Field(void) const {
	if (!fField) ExceptionT::GeneralFail("FBC_ControllerT::Field", "pointer not set");
	return *fField;
}

/* return the FieldSupportT or throw exception if not set */
inline const FieldSupportT& FBC_ControllerT::FieldSupport(void) const {
	if (!fFieldSupport) ExceptionT::GeneralFail("FBC_ControllerT::FieldSupport", "pointer not set");
	return *fFieldSupport;
}

} /* namespace Tahoe */

#endif /* _FBC_CONTROLLER_T_H_ */
