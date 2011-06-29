/* $Id: BTreeNodeT.h,v 1.2 2002-11-18 10:00:21 paklein Exp $ */
#ifndef _B_TREE_NODE_T_H_
#define _B_TREE_NODE_T_H_

namespace Tahoe {

template <class TYPE> class BinaryTreeT;

/** container class for BinaryTreeT. 
 * \note the TYPE stored in the list should have an appropriate
 * copy constructor. */
template <class TYPE>
class BTreeNodeT
{
	friend class BinaryTreeT<TYPE>;

public:

	/** constructor */
	BTreeNodeT(const TYPE &value);
	 	 	
	/** return data from the node */
	const TYPE& Value(void) const;

	/** pointer to the left node */
	BTreeNodeT<TYPE>* Left(void) const { return fLeft; };

	/** pointer to the right node */
	BTreeNodeT<TYPE>* Right(void) const { return fRight; };
	 	 	
private:

	TYPE fValue;
	BTreeNodeT*	fLeft;
	BTreeNodeT*	fRight;
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructor */
template <class TYPE>
inline BTreeNodeT<TYPE>::BTreeNodeT(const TYPE &value):
	fValue(value),
	fLeft(NULL),
	fRight(NULL)	
{
	
}

/* return data from the node */
template <class TYPE>
inline const TYPE& BTreeNodeT<TYPE>::Value(void) const
{
	return fValue;
}

} // namespace Tahoe 
#endif /* _B_TREE_NODE_T_H_ */
