/* $Id: UnConnectedRodT.h,v 1.7 2004-01-05 07:19:56 paklein Exp $ */
/* created: paklein (04/05/1997) */

#ifndef _UNCONN_ROD_T_H_
#define _UNCONN_ROD_T_H_

/* base class */
#include "RodT.h"

namespace Tahoe {

/** Interface for a rod element group that connects itself based on the
 * nodes placed in the group. All the rods in the group are assumed to
 * be identical, ie. only 1 material set may be specified in the input. */
class UnConnectedRodT: public RodT
{
public:

	/* constructor */
	UnConnectedRodT(const ElementSupportT& support, const FieldT& field);

	/* apply pre-conditions at the current time step.  Signal
	 * all listeners that the time has just been incremented */
	virtual void InitStep(void);

	/* resets to the last converged solution */
	virtual GlobalT::RelaxCodeT ResetStep(void);

	/* element level reconfiguration for the current solution */
	virtual GlobalT::RelaxCodeT RelaxSystem(void);
			
protected: /* for derived classes only */

	/* print element group data */
	virtual void PrintControlData(ostream& out) const;
	 			
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/** return true if connectivities are changing */
	virtual bool ChangingGeometry(void) const { return true; };

private:

	/* configure element list data */
	void ConfigureElementData(void);

	/* print data */
	void PrintConnectivityData(ostream& out);

protected:

	/** changing connectivities */
	//iArray2DT fRodConnectivities;

private:

	/* neighbor resetting increment */
	int fReconnectInc;

	/* neighbor searching parameters */
	int     fMaxNeighborCount;
	double	fNeighborDist;
	int		fNumNodesUsed; // equals -1 for all nodes in the NodeManagerT
	//need to store NodesUsed for reconnection. Not implemented
	
	/* runtime data */
	int fReconnectCount;
};

} // namespace Tahoe 
#endif /* _UNCONN_ROD_T_H_ */
