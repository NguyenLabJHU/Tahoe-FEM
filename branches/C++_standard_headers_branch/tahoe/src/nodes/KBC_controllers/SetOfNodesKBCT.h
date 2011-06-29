/* $Id: SetOfNodesKBCT.h,v 1.3 2004-07-15 08:31:21 paklein Exp $ */
#ifndef _SET_OF_NODES_KBC_T_H_
#define _SET_OF_NODES_KBC_T_H_

/* base class */
#include "KBC_ControllerT.h"

/* direct members */
#include "ScheduleT.h"
#include "iArrayT.h"
#include "BasicFieldT.h"

namespace Tahoe {

/** Nodes whose velocities are scaled to have zero net momentum
 * and kinetic energies that sum to 3/2 N k_B T
 */
class SetOfNodesKBCT: public KBC_ControllerT
{
public:	

	/** constructor */
	SetOfNodesKBCT(const BasicSupportT& support, BasicFieldT& field);

	/** destructor */
 	~SetOfNodesKBCT(void);

	/** initialize data. Must be called immediately after construction */
	virtual void Initialize(ifstreamT& in);

	virtual void InitialCondition(void) { } ;

protected:
	
	void SetBCCards(void);
	/*@}*/

protected:

	/** data for schedule for the KBC */
	const ScheduleT* fSchedule;
	int fScheduleNum;
	double fScale;
	
	/** nodes in spatial regions that are being controlled */
	iArrayT fNodes;
	
	/** re-check for nodes in region every this many timesteps */
	int nIncs;

};

} // namespace Tahoe 
#endif /* _SET_OF_NODES_KBC_T_H_ */
