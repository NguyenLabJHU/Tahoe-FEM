/* $Id: VirtualSWDC.h,v 1.3 2002-07-02 19:55:28 cjkimme Exp $ */
/* created: paklein (05/05/1997) */
#ifndef _VIRTUAL_SWDC_H_
#define _VIRTUAL_SWDC_H_

/* base class */
#include "SWDiamondT.h"


namespace Tahoe {

class VirtualSWDC: public SWDiamondT
{
public:

	/* constructor */
	VirtualSWDC(const ElementSupportT& support, const FieldT& field);

	/* append element equations numbers to the list */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/* appends group connectivities to the array */
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;
		 			  	
protected: /* for derived classes only */
	 			
	/* element data */
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

private:

	/* blind substitution of all virtual node pairs */
	void SwapVirtualNodes(iArray2DT& elnodelist) ; //not quite
	void SwapVirtualNodes2(iArray2DT& elnodelist); //even less quite

private:

	iArray2DT fPeriodicNodes_3Body;
	iArray2DT fVNodePairs;
};

} // namespace Tahoe 
#endif /* _VIRTUAL_SWDC_H_ */
