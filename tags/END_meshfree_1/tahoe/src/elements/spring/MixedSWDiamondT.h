/* $Id: MixedSWDiamondT.h,v 1.1.1.1 2001-01-29 08:20:38 paklein Exp $ */
/* created: paklein (03/22/1997)                                          */
/* Interface for heterogeneous diamond cubic lattice                      */

#ifndef _MIXED_SWDIAMOND_T_H_
#define _MIXED_SWDIAMOND_T_H_

/* base class */
#include "SWDiamondT.h"

/* direct members */
#include "SWDataT.h"

/* forward declarations */
class LoadTime;

class MixedSWDiamondT: public SWDiamondT
{
public:

	/*
	 * constructor
	 */
	MixedSWDiamondT(FEManagerT& fe_manager);

	/*
	 * Apply pre-conditions at the current time step.  Signal
	 * all listeners that the time has just been incremented.
	 */
	virtual void InitStep(void);

protected:

	/*
	 * Print element group data.
	 */
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
	int 		fLTfNum;
	LoadTime* 	fLTfPtr;

	/* material set list */
	ArrayT<SWDataT> fSWDataList;
	int		        fCurrMatType;

	/* material tags for every node */
	iArrayT fNodeTypes;	// (1...numnodes) 			  	
};

#endif /* _MIXED_SWDIAMOND_T_H_ */
