/* $Id: ElementBlockDataT.h,v 1.4 2003-08-14 05:32:30 paklein Exp $ */
#ifndef _ELEM_BLOCK_DATA_T_H_
#define _ELEM_BLOCK_DATA_T_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "StringT.h"

namespace Tahoe {

/** container for element block information. Used by ElementBaseT to
 * store information about multiple element blocks within an element
 * group. */
class ElementBlockDataT: public ParameterInterfaceT
{
public:
	
	/** constructor */
	ElementBlockDataT(void);

	/** copy constructor */
	ElementBlockDataT(ElementBlockDataT& source);

	/** set all field */
	void Set(const StringT& ID, int start, int dim, int material);
	
	/** set the number of elements in the block */
	void SetDimension(int dim) { fDimension = dim; };

	/** element block ID */
	const StringT& ID(void) const { return fID; };
	
	/** number of block's first element with the group */
	int StartNumber(void) const { return fStartNum; };
	
	/** number of elements in the block */
	int Dimension(void) const { return fDimension; };
	
	/** material ID */
	int MaterialID(void) const { return fMaterial; };
	
	/** assigmnent operator */
	ElementBlockDataT& operator=(const ElementBlockDataT& rhs);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;
	/*@}*/
	
private:

	/** block ID */
	StringT fID;
	
	/** start number within the element group */
	int fStartNum;
	
	/** number of elements in the block */
	int fDimension;
	
	/** material identifier */
	int fMaterial;
};

/* inlines */
inline ElementBlockDataT& ElementBlockDataT::operator=(const ElementBlockDataT& rhs)
{
	Set(rhs.fID, rhs.fStartNum, rhs.fDimension, rhs.fMaterial);
	return *this;
}

} // namespace Tahoe 
#endif /* _ELEM_BLOCK_DATA_T_H_ */
