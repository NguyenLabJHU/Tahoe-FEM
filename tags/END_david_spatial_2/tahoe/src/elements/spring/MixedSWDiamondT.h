/* $Id: MixedSWDiamondT.h,v 1.4.86.1 2005-05-27 19:55:18 paklein Exp $ */
/* created: paklein (03/22/1997) */

#ifndef _MIXED_SWDIAMOND_T_H_
#define _MIXED_SWDIAMOND_T_H_

/* base class */
#include "SWDiamondT.h"

/* direct members */
#include "SWDataT.h"

namespace Tahoe {

/* forward declarations */
class ScheduleT;

class MixedSWDiamondT: public SWDiamondT
{
public:

	/** constructor */
	MixedSWDiamondT(const ElementSupportT& support, const FieldT& field);

	/** init new time increment */
	virtual GlobalT::InitStatusT InitStep(void);

protected:

	/** print element group data */
	virtual void PrintControlData(ostream& out) const;
	 			
	/* element data */
	virtual void ReadMaterialData(ifstreamT& in);	
	virtual void WriteMaterialData(ostream& out) const;
	virtual void EchoConnectivityData(ifstreamT& in, ostream& out);

	/* element list increment */
	virtual bool Next2Body(void);
	virtual bool Next3Body(void);

private:

	/* echo node type tags */
	void EchoNodeTags(istream& in, ostream& out);

	/* copy material properties from the specified set */
	void CopyMaterialData(int setnum);
	void Mix2Body(int m1, int m2);
	void Mix3Body(int m1, int m_mid, int m2);
		 			
private:
	
	/* variation LTf */
	int fLTfNum;
	const ScheduleT* fLTfPtr;

	/* material set list */
	ArrayT<SWDataT> fSWDataList;
	int		        fCurrMatType;

	/* material tags for every node */
	iArrayT fNodeTypes;	// (1...numnodes) 			  	
};

} // namespace Tahoe 
#endif /* _MIXED_SWDIAMOND_T_H_ */
