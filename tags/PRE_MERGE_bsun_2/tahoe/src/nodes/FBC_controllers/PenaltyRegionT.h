/* $Id: PenaltyRegionT.h,v 1.6 2003-10-04 19:14:05 paklein Exp $ */
/* created: paklein (04/30/1998) */

#ifndef _PENALTY_REGION_T_H_
#define _PENALTY_REGION_T_H_

/* base class */
#include "FBC_ControllerT.h"

/* direct members */
#include "iArrayT.h"
#include "dArrayT.h"
#include "dArray2DT.h"
#include "iArray2DT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

/** base class for moving rigid, penalty regions. contact nodes
 * that enter the region are expelled by a quadratic penetration
 * potential. derived classes are responsilble for computing
 * the penetration depth and reaction for based on the geometry
 * of the region. */
class PenaltyRegionT: public FBC_ControllerT
{
public:

	/** motion control codes */
	enum MotionCodeT {
		kConstantVelocity = 0, /**< on change in velocity */
		         kImpulse = 1, /**< region slows with contact impulse */
		        kSchedule = 2  /**< velocity follows schedule function */
			};

	/* constructor */
	PenaltyRegionT(FEManagerT& fe_manager, int group, const iArray2DT& eqnos,
		const dArray2DT& coords, const dArray2DT& displ, const dArray2DT* vels);

	/* input processing */
	virtual void EchoData(ifstreamT& in, ostream& out);

	/* initialize data */
	virtual void Initialize(void);

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

	/* initialize/finalize step */
	virtual void InitStep(void);
	virtual void CloseStep(void);

	/* reset displacements (and configuration to the last known solution) */
	virtual void Reset(void);

	/** \name writing results */
	/*@{*/
	/** register data for output */
	virtual void RegisterOutput(void);

	/** write results */
	virtual void WriteOutput(ostream& out) const;
	/*@}*/

private:

	/* accumulate the contact force vector fContactForce */
	virtual void ComputeContactForce(double kforce) = 0;

protected:

	/** \name references to NodeManagerT data */
	/*@{*/
	const iArray2DT& rEqnos;  /**< nodal equation numbers */
	const dArray2DT& rCoords; /**< nodal coordinates */
	const dArray2DT& rDisp;   /**< nodal displacement */
	const dArray2DT* pVels;   /**< nodal velocities */
	/*@}*/

	/** \name wall input parameters */
	/*@{*/
	dArrayT fx0;             /**< initial position */
	dArrayT fv0;             /**< initial velocity */
	double fk;               /**< penalty stiffness */
	int	   fSlow;            /**< 1 if the region slows from collisions */
	double fMass;            /**< mass of the region */
	const ScheduleT* fLTf;   /**< NULL if there is no time dependence */
	int    fNumContactNodes; /**< number of contact nodes */
	/*@}*/

	/** \name state variables */
	/*@{*/
	dArrayT fx;     /**< position */
	dArrayT fv;     /**< velocity */
	dArrayT fxlast; /**< last converged position */
	dArrayT fvlast; /**< last converged velocity */
	/*@}*/

	/** \name contact force node and equation numbers */
	/*@{*/
	iArrayT fContactNodes;
	iArrayT fContactEqnos;

	/** shallow version of PenaltyRegionT::fContactForce2D */
	dArrayT fContactForce;

	/** array of signed gaps where gap < 0.0 implies contact */
	dArrayT fGap;

	/** dArray2DT copy of the force */
	dArray2DT fContactForce2D;
	/*@}*/

	/* workspace */
	dArrayT fTempNumNodes; // temp space length = fNumContactNodes

	/** \name writing results */
	/*@{*/	
	/** output ID */
	int fOutputID;
	
	/** "connectivities" for output, just alias of PenaltyRegionT::fContactNodes */
	iArray2DT fContactNodes2D;
	/*@}*/	
};

} // namespace Tahoe 
#endif /* _PENALTY_REGION_T_H_ */
