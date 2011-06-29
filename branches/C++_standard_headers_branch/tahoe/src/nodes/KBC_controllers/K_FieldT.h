/* $Id: K_FieldT.h,v 1.12 2009-05-21 22:30:27 tdnguye Exp $ */
/* created: paklein (09/05/2000) */
#ifndef _K_FIELD_T_H_
#define _K_FIELD_T_H_

#include "ElementsConfig.h"
#ifdef CONTINUUM_ELEMENT

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "dArrayT.h"
#include "iArrayT.h"
#include "ScheduleT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class ElementBaseT;
class IsotropicT;
class SolidMaterialT;

/** K-field displacements */
class K_FieldT: public KBC_ControllerT
{
public:

	/** tip tracking methods. Define how the crack tip is determined from the
	 * nodal values returned by the crack tip tracking group set by K_FieldT::fNearTipGroupNum */
	enum TrackingCodeT {
		 kMaximum = 0, /**< location of the maximum value */
	   kThreshold = 1  /**< location of the first value exceeding a threshold */
	};

	/* constructor */
	K_FieldT(const BasicSupportT& support);

	/* initial condition/restart functions
	 *
	 * Set to initial conditions.  The restart functions
	 * should read/write any data that overrides the default
	 * values */
	virtual void InitialCondition(void);
	virtual void ReadRestart(ifstreamT& in);
	virtual void WriteRestart(ofstreamT& out) const;

	/* initialize/finalize/reset step */
	virtual void InitStep(void);
	virtual void CloseStep(void);
	virtual void Reset(void);

	/* returns true if the internal force has been changed since
	 * the last time step */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);

	/* output current configuration */
	virtual void WriteOutput(ostream& out) const;

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** extract elastic constants */
	void ResolveElasticProperties(const ParameterListT& list, int& group_number, int& material_number, double& mu, double& nu, double& kappa) const;

	/* determine the new tip coordinates */
	void GetNewTipCoordinates(dArrayT& tip_coords);

	/* resolve element info to isotropic material */
	void ResolveMaterialReference(int element_group, int material_num,
		const IsotropicT** piso, const SolidMaterialT** pmat) const;

	/** compute K-field displacement factors. Recompute the asymptotic displacement
	 * field as a function of the current values of K_FieldT::fmu and K_FieldT::fkappa. */
	virtual void ComputeDisplacementFactors(const dArrayT& tip_coords);
	
	/* set BC cards with current displacement field */
	virtual void SetBCCards(void);
	
protected:

	/** \name K-field specifications */
	/*@{*/
	double fK1;
	double fK2;
	/*@}*/

	/** crack tip coordinates */
	dArrayT fInitTipCoords;
	
	/** crack extension parameters */
	dArrayT fGrowthDirection;

	/** \name crack tip tracking parameters */
	/*@{*/
	/** near tip group or -1 to disable any tracking */
	int fNearTipGroupNum;

	/** near tip output variable */
	StringT fNearTipOutputVariable;
	
	/** nodal output code from tip group used to locate crack tip */
	int fNearTipOutputCode;

	/** value within the output variables to locate tip */
	int fTipColumnNum;

	/** tip tracking method */
	TrackingCodeT fTrackingCode;

	/** data used for tip tracking */
	dArrayT fTrackingParameters;
	/*@}*/

	/** \name crack extension limiters */
	/*@{*/
	/** total extension during a single time increment */
	double fMaxGrowthDistance;

	/** maximum number of relaxation steps within a single time increment */
	int fMaxGrowthSteps;
	/*@}*/
		
	/* BC nodes */
	ArrayT<StringT> fID_List;
	iArrayT fNodes;
	
	/** \name elastic constants */
	/*@{*/
	double fmu; /**< shear modulus */
	double fnu; /**< Poisson's ratio */
	double fkappa; /**< function of nu */
	
	int fGroupNumber;
	int fMaterialNumber;
	/*@}*/

	/* runtime data */
	ScheduleT fDummySchedule;
	const ScheduleT* fLTf1;
	const ScheduleT* fLTf2;   	
	dArray2DT fK1Disp;
	dArray2DT fK2Disp;
	int fGrowthCount;

	/* tip coordinates */
	dArrayT fTipCoords;
	dArrayT fLastTipCoords;
};

} /* namespace Tahoe */

#endif /* CONTINUUM_ELEMENT */
#endif /* _K_FIELD_T_H_ */
