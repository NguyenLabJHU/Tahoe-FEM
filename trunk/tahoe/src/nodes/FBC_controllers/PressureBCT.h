/* $Id: PressureBCT.h,v 1.2 2006-11-16 00:59:52 r-jones Exp $ */

#ifndef _PRESSURE_BC_T_H_
#define _PRESSURE_BC_T_H_

/* base class */
#include "FBC_ControllerT.h"
#if 0
#include "DOFElementT.h"
#endif

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
class XDOF_ManagerT;


#if 0
class PressureBCT: public FBC_ControllerT, public DOFElementT
#endif
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
  virtual void Connectivities(AutoArrayT<const iArray2DT*>& connects_1,
    AutoArrayT<const RaggedArray2DT<int>*>& connects_2,
    AutoArrayT<const iArray2DT*>& equivalent_nodes) const;


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

	enum ControlTypesT { kPressureControl = 0,
	                     kVolumeControl = 1 };

#if 0
  /** \name implementation of the DOFElementT interface */
  /*@{*/
  /* returns the array for the DOF tags needed for the current config */
  virtual void SetDOFTags(void);
  virtual iArrayT& DOFTags(int tag_set);

  /* generate nodal connectivities - does nothing here */
  virtual void GenerateElementData(void);

  /* return the contact elements */
  virtual const iArray2DT& DOFConnects(int tag_set) const;

  /* restore the DOF values to the last converged solution */
  virtual void ResetDOF(dArray2DT& DOF, int tag_set) const;

  /* returns 1 if group needs to reconfigure DOF's, else 0 */
  virtual int Reconfigure(void);

	/** restore any state data to the previous converged state */
	virtual void ResetState(void) { };

	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const;
	/*@}*/
#endif

private:
	void InputSideSets(const ParameterListT& list, 
			GeometryT::CodeT geom_code, 
			iArray2DT& facets);
	void ComputeVolume(DomainIntegrationT& domain, dArray2DT& coord,
			                double& volume, double& area);
	void ComputeForce(DomainIntegrationT& domain, dArray2DT& coord,
			                dArray2DT& force);
	void ComputeStiffness(DomainIntegrationT& domain, dArray2DT& coord, 
			                ElementMatrixT& stiffness);
	void ComputeVolumeStiffness(DomainIntegrationT& domain, dArray2DT& coord, 
			                dArray2DT& delV);
	int Connectivities(iArray2DT& conn);


protected:

	ArrayT<iArray2DT> fSurfaces;
	ArrayT<DomainIntegrationT*> fDomains;
	const ScheduleT* fSchedule ; /**< NULL if there is no time dependence */
	double fScheduleScale; /**< schedule scaling to get pressure (or change in volumee */
	int fNumNodes; /**< number of nodes */
	int fNumSurfaces; /**< number of surfaces */
	int fOutputID; /** output ID */
	int fControlType; /** control pressure or volume */
	double fPenalty; /** penalty for volume control */
	int fndir; /** gross normal direction for surfae, used to compute volume */
	int fnsd; /** number of spatial dimensions */
	double fPressure; /** pressure on surface */
	double fVolume0; /** initial enclosed volume in normal direction */
	double fVolume; /** enclosed volume in normal direction */
	double fArea; /** projected area in normal direction */
	dArrayT fReaction; /** reaction forces */
	iArray2DT fConnectivities;
	iArray2DT fEquationNumbers;
	iArrayT fMultiplierTags;
	iArray2DT fMultiplierConnects;
	
};

} // namespace Tahoe 
#endif /* _PRESSURE_BC_T_H_ */
