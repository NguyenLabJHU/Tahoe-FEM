/* $Id: ConstantStressNodesT.h,v 1.1.2.1 2003-09-18 21:03:44 cjkimme Exp $ */
#ifndef _CONSTANT_STRESS_NODES_T_H_
#define _COSNTANT_STRESS_NODES_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "ScheduleT.h"
#include "iArrayT.h"
#include "BasicFieldT.h"

namespace Tahoe {

/* forward declarations */
class BasicFieldT;
class RandomNumberT;

/** Nodes whose velocities are scaled to have zero net momentum
 * and kinetic energies that sum to 3/2 N k_B T
 */
class ConstantStressNodesT: public KBC_ControllerT
{
public:	

	/** constructor */
	ConstantStressNodesT(NodeManagerT& node_manager, BasicFieldT& field);

	/** destructor */
 	~ConstantStressNodesT(void);

	/** initialize data. Must be called immediately after construction */
	virtual void Initialize(ifstreamT& in);

	/** write class parameters */
	void WriteParameters(ostream& out) const;

	/** do at start of timestep */
	virtual void InitStep(void);

	/** Initialize to appropriate temperature */
	virtual void InitialCondition(void);
	
	virtual bool IsICController(void) { return true; }

protected:
	
	void SetBCCards(void);
	
	/** apply the update to the solution. Corrects velocities. */
	virtual void Update(const dArrayT& update);
	/*@}*/

protected:

	/** the field */
	BasicFieldT& fField;

	/** True if controller only used for IC */
	bool qIConly;
	
	/** flag to let this controller only influence ICs */
	bool qFirstTime;
	
	/** true if allNodes needs to be initialized or rescaled */
	bool qAllNodes;

	/** adjust pressure every fIncs timesteps */
	int fIncs, fIncCt;

	/** evolution controlled by a schedule */
	const ScheduleT* fTempSchedule;
	int fnumTempSchedule;
	double fTempScale;
	
	/** temperature schedule is not the BC value. Need a dummy schedule, too */
	ScheduleT fDummySchedule;
	
	/** constrained nodes */
	/*@{*/
	/** id list for the \e node sets */
	ArrayT<StringT> fNodeIds;
	iArrayT fNodes;

	/** assuming all nodes have same mass */
	double fMass;
	
	/** initial pressure */
	double fP_0;
};

} // namespace Tahoe 
#endif /* _SCALED_VELOCITY_NODES_T_H_ */
