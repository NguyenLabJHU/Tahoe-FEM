/* $Id: VirtualRodT.h,v 1.1.1.1 2001-01-29 08:20:34 paklein Exp $ */
/* created: paklein (05/01/1997)                                          */
/* UnConnectedRodT plus virtual elements for periodic boundary            */
/* conditions.                                                            */

#ifndef _UNCON_VROD_T_H_
#define _UNCON_VROD_T_H_

/* base class */
#include "UnConnectedRodT.h"

class VirtualRodT: public UnConnectedRodT
{
public:

	/* constructor */
	VirtualRodT(FEManagerT& fe_manager);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

protected:

	/* element data */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

private:

	/* shuffle node numbers */
	void SwapVirtualNodes(iArray2DT& elnodelist) const;
	void SwapVirtualNodes2(iArray2DT& elnodelist); //blind swap of all virtual node pairs

private:

	iArray2DT	fVNodeTriplets; //needs to be triplet in order to isolate
	                            //those 2 body interactions to which the
	                            //virtual node pair applies
};

#endif /* _UNCON_VROD_T_H_ */
