/* $Id: PressureBCT.h,v 1.1 2006-06-19 15:25:34 r-jones Exp $ */

#ifndef _PRESSURE_BC_T_H_
#define _PRESSURE_BC_T_H_

/* base class */
#include "FBC_ControllerT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dMatrixT.h"
#include "DomainIntegrationT.h"
#include "ElementMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

class PressureBCT: public FBC_ControllerT
{
public:

	/** constructor */
	PressureBCT(void);

	/** destructor */
	virtual ~PressureBCT(void);

	/* form of tangent matrix */
	virtual GlobalT::SystemTypeT TangentType(void) const;

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/* initial condition/restart functions
	 *
	 * Set to initial conditions.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	virtual void InitialCondition(void);
	virtual void ReadRestart(istream& in);
	virtual void WriteRestart(ostream& out) const;

	/* apply force */
	virtual void ApplyRHS(void);

	/* tangent */
	virtual void ApplyLHS(GlobalT::SystemTypeT sys_type);

	/* initialize/finalize step */
	virtual void InitStep(void);
	virtual void CloseStep(void);

	/* reset displacements (and configuration to the last known solution) */
	virtual void Reset(void);

	/** returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/** \name writing results */
	/*@{*/
	/** register data for output */
	virtual void RegisterOutput(void);

	/** write results */
	virtual void WriteOutput(ostream& out) const;
	/*@}*/

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

private:
	void InputSideSets(const ParameterListT& list, 
			GeometryT::CodeT geom_code, 
			iArray2DT& facets);
	void Compute_Force(DomainIntegrationT& domain, dArray2DT& coord,
			                dArray2DT& force);
	void Compute_Stiffness(DomainIntegrationT& domain, dArray2DT& coord, 
			                ElementMatrixT& stiffness);


protected:

        ArrayT<iArray2DT> fSurfaces;
        ArrayT<DomainIntegrationT*> fDomains;
	const ScheduleT* fSchedule ; /**< NULL if there is no time dependence */
	double fpscale; /**< pressure scale */
	int fNumNodes; /**< number of nodes */
	int fNumSurfaces; /**< number of surfaces */
	int fOutputID; /** output ID */
	
};

} // namespace Tahoe 
#endif /* _PRESSURE_BC_T_H_ */
