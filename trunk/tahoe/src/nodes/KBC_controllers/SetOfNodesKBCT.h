/* $Id: SetOfNodesKBCT.h,v 1.1 2003-05-28 17:44:52 cjkimme Exp $ */
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
	SetOfNodesKBCT(NodeManagerT& node_manager, BasicFieldT& field);

	/** destructor */
 	~SetOfNodesKBCT(void);

	/** initialize data. Must be called immediately after construction */
	virtual void Initialize(ifstreamT& in);

	virtual void InitialCondition(void) { } ;
	
	virtual void WriteParameters(ostream& out) const;

protected:
	
	void SetBCCards(void);
	/*@}*/

protected:

	/** temperature evolution controlled by a schedule */
	const ScheduleT* fSchedule;
	int fScheduleNum;
	double fScale;

};

} // namespace Tahoe 
#endif /* _SET_OF_NODES_KBC_T_H_ */
