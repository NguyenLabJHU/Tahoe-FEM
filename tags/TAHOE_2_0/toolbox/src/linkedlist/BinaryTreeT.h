/* $Id: BinaryTreeT.h,v 1.4 2003-03-09 20:37:19 paklein Exp $ */
#ifndef _BINARY_TREE_T_H_
#define _BINARY_TREE_T_H_

#include "ExceptionT.h"

/* direct members */
#include "BTreeNodeT.h"
#include "ArrayT.h"

namespace Tahoe {

/** basic templated binary tree class
 * \note the TYPE stored in the list should have an appropriate
 * copy constructor, an assignment "=", and inequality operators "<" and ">" */
template <class TYPE>
class BinaryTreeT
{
public:

	/** constructor */
	BinaryTreeT(void);

	/** destructor */
	~BinaryTreeT(void);

	/** override or reset the function used to compare tree values during calls
	 * to BinaryTreeT::Find. It should return
	 * <ul> 
	 * <li> <0 : tree_value < test_value
	 * <li> >0 : tree_value > test_value
	 * <li>  0 : tree_value == test_value
	 * </ul>
	 * Note that the order of the arguments matters in the comparison. Pass NULL to
	 * clear the comparison function and use the default comparison operators < and >.
	 * \param compare the comparison function
	 * \param tree_value value in currently in the tree
	 * \param test_value value being compared to the tree
	 */
	void SetCompareFunction(int (*compare)(const TYPE& tree_value, const TYPE& test_value));

	/** \name insert values into the tree */
	/*@{*/
	/** insert single value */
	void Insert(const TYPE& value);

	/** insert an array of values */
	void Insert(const ArrayT<TYPE>& values);

	/** insert a tree of values */
	void Insert(const BinaryTreeT<TYPE>& values);
	/*@}*/

	/** insert the values if unique */
	/*@{*/
	/** insert single value. Return true if the value is inserted. */
	bool InsertUnique(const TYPE& value);

	/** insert an array of values. Return number of values inserted. */
	int InsertUnique(const ArrayT<TYPE>& values);

	/** insert a tree of values */
	int InsertUnique(const BinaryTreeT<TYPE>& values);
	/*@}*/
				
	/** clears tree contents */
	void Clear(void);	
		
	/** returns the number of values in the list */	
	int Size(void) const { return fSize; };
	
	/** return true if the tree contains the given value */
	bool HasValue(const TYPE& value) const { return Find(value) != NULL; };
	
	/** \name return as array */
	/*@{*/
	/** return the tree values in ascending order. The destination array is
	 * dimensioned during the call. */
	void Ascending(ArrayT<TYPE>& array) const;

	/** return the tree values in descending order. The destination array is
	 * dimensioned during the call. */
	void Descending(ArrayT<TYPE>& array) const;
	/*@}*/
	
	/** return the tree node with the given value or NULL if not present */
	BTreeNodeT<TYPE>* Find(const TYPE& value) const;

private:

	/** \name not allowed */
	/*@{*/
	/** copy constructor */
	BinaryTreeT(const BinaryTreeT&);
		
	/** assignment operator */
	BinaryTreeT& operator=(const BinaryTreeT&);
	/*@}*/

	/** \name recursive tree operators */
	/*@{*/
	/** recursive function to perform the clear */
	void DoClear(BTreeNodeT<TYPE>* node);

	/** insert the node in the tree
	 * \param node node in this tree
	 * \param value value to insert */
	void DoInsert(BTreeNodeT<TYPE>* node, const TYPE& value);

	/** insert the source tree in the tree
	 * \param node node in the source tree */
	void DoInsert(BTreeNodeT<TYPE>* node);

	/** insert the node in the tree. Return true if inserted
	 * \param node node in this tree
	 * \param value value to insert */	
	bool DoInsertUnique(BTreeNodeT<TYPE>* node, const TYPE& value);

	/** insert unique values from the source tree
	 * \param node node in the source tree */
	int DoInsertUnique(BTreeNodeT<TYPE>* node);

	/** return true if a matching value is found */
	BTreeNodeT<TYPE>* DoFind(BTreeNodeT<TYPE>* node, const TYPE& value) const;

	/** gather values in ascending order */
	void DoAscending(BTreeNodeT<TYPE>* node, ArrayT<TYPE>& array, int& dex) const;

	/** gaher values in descending order */
	void DoDescending(BTreeNodeT<TYPE>* node, ArrayT<TYPE>& array, int& dex) const;
	/*@}*/

	/** \name creating and destroying nodes.
	 * Keep BinaryTreeT::fSize up to date */
	/*@{*/
	/** return a pointer to a new node containing the given value */
	BTreeNodeT<TYPE>* NewNode(const TYPE& value);

	/** delete the given node. Does nothing if node is NULL. */
	void FreeNode(BTreeNodeT<TYPE>* node);
	/*@}*/

private:

	/** number of nodes in the tree */
	int fSize;

	/** pointer to the first node in the list */
	BTreeNodeT<TYPE>* fRoot;
	
	/** comparison function used by BinaryTree::Find */
	int (*fCompare)(const TYPE&, const TYPE&);
};

/*************************************************************************
* Implementation
*************************************************************************/

template <class TYPE>
BinaryTreeT<TYPE>::BinaryTreeT(void):
	fSize(0),
	fRoot(NULL),
	fCompare(NULL)
{

}
	
template <class TYPE>
BinaryTreeT<TYPE>::~BinaryTreeT(void) { Clear(); }

template <class TYPE>
void BinaryTreeT<TYPE>::SetCompareFunction(int (*compare)(const TYPE& tree_value, const TYPE& test_value))
{
	fCompare = compare;
}

/* insert single value */
template <class TYPE>
inline void BinaryTreeT<TYPE>::Insert(const TYPE& value)
{
	/* empty */
	if (!fRoot)
		fRoot = NewNode(value);
	else 
		DoInsert(fRoot, value);
}

/* insert an array of values */
template <class TYPE>
inline void BinaryTreeT<TYPE>::Insert(const ArrayT<TYPE>& values)
{
	for (int i = 0; i < values.Length(); i++)
		Insert(values[i]);
}

/* insert a tree of values */
template <class TYPE>
void BinaryTreeT<TYPE>::Insert(const BinaryTreeT<TYPE>& values)
{
	/* start running through the tree */
	if (values.fRoot) DoInsert(values.fRoot);
}

/* insert single value */
template <class TYPE>
inline bool BinaryTreeT<TYPE>::InsertUnique(const TYPE& value)
{
	/* empty */
	if (!fRoot)
	{
		fRoot = NewNode(value);
		return true;
	}
	else 
		return DoInsertUnique(fRoot, value);
}

/* insert an array of values */
template <class TYPE>
inline int BinaryTreeT<TYPE>::InsertUnique(const ArrayT<TYPE>& values)
{
	int count = 0;
	for (int i = 0; i < values.Length(); i++)
		if (InsertUnique(values[i]))
			count++;
	return count;
}

/* insert a tree of values */
template <class TYPE>
inline int BinaryTreeT<TYPE>::InsertUnique(const BinaryTreeT<TYPE>& values)
{
	/* skip this */
	if (&values == this)
		return 0;
	/* start running through the tree */
	else if (values.fRoot) 
		return DoInsertUnique(values.fRoot);
	else
		return 0;
}

/* free all memory */
template <class TYPE>
inline void BinaryTreeT<TYPE>::Clear(void)
{ 
	DoClear(fRoot);
	fRoot = NULL;
}

/* return tree in ascending order */
template <class TYPE>
void BinaryTreeT<TYPE>::Ascending(ArrayT<TYPE>& array) const
{
	array.Dimension(Size());
	if (Size() > 0)
	{
		int dex = 0;
		DoAscending(fRoot, array, dex);
		if (dex != array.Length()) ExceptionT::GeneralFail();
	}
}

/* return the tree values in descending order */
template <class TYPE>
void BinaryTreeT<TYPE>::Descending(ArrayT<TYPE>& array) const
{
	array.Dimension(Size());
	if (Size() > 0)
	{
		int dex = 0;
		DoDescending(fRoot, array, dex);
		if (dex != array.Length()) ExceptionT::GeneralFail();
	}
}

/* return true if the tree contains the given value */
template <class TYPE>
BTreeNodeT<TYPE>* BinaryTreeT<TYPE>::Find(const TYPE& value) const
{
	if (!fRoot)
		return NULL;
	else
		return DoFind(fRoot, value);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* free all memory */
template <class TYPE>
void BinaryTreeT<TYPE>::DoClear(BTreeNodeT<TYPE>* node) 
{
	if (node != NULL)
	{
		DoClear(node->fLeft); 
		DoClear(node->fRight);
		FreeNode(node);
	}
}

/* insert the node in the tree */
template <class TYPE>
void BinaryTreeT<TYPE>::DoInsert(BTreeNodeT<TYPE>* node, const TYPE& value)
{
	/* look left */
	if (value < node->fValue)
	{
		if (node->fLeft)
			DoInsert(node->fLeft, value);
		else
			node->fLeft = NewNode(value);
	}
	/* look right */
	else
	{
		if (node->fRight)
			DoInsert(node->fRight, value);
		else
			node->fRight = NewNode(value);
	}
}

/* insert the node in the tree */
template <class TYPE>
void BinaryTreeT<TYPE>::DoInsert(BTreeNodeT<TYPE>* node)
{
	/* collect lower tree */
	if (node->fLeft) DoInsert(node->fLeft);

	/* collect this */
	Insert(node->fValue);
	
	/* collect upper tree */
	if (node->fRight) DoInsert(node->fRight);
}

/* insert if unique */
template <class TYPE>
bool BinaryTreeT<TYPE>::DoInsertUnique(BTreeNodeT<TYPE>* node, const TYPE& value)
{
	/* look left */
	if (value < node->fValue)
	{
		if (node->fLeft)
			return DoInsertUnique(node->fLeft, value);
		else
		{
			node->fLeft = NewNode(value);
			return true;
		}
	}
	/* look right */
	else if (value > node->fValue)
	{
		if (node->fRight)
			return DoInsertUnique(node->fRight, value);
		else
		{
			node->fRight = NewNode(value);
			return true;
		}
	}
	/* same as node->fValue */
	return false;
}

/* insert the node in the tree */
template <class TYPE>
int BinaryTreeT<TYPE>::DoInsertUnique(BTreeNodeT<TYPE>* node)
{
	int count = 0;

	/* collect lower tree */
	if (node->fLeft) count += DoInsertUnique(node->fLeft);

	/* collect this */
	if (InsertUnique(node->fValue)) count++;
	
	/* collect upper tree */
	if (node->fRight) count += DoInsertUnique(node->fRight);
	
	return count;
}

/* return true if found */
template <class TYPE>
BTreeNodeT<TYPE>* BinaryTreeT<TYPE>::DoFind(BTreeNodeT<TYPE>* node, const TYPE& value) const
{
	/* user-defined test function */
	if (fCompare)
	{
		int compare = fCompare(node->fValue, value);

		/* look left */
		if (compare > 0)
		{
			if (node->fLeft)
				return DoFind(node->fLeft, value);
			else
				return NULL;
		}
		/* look right */
		if (compare < 0)
		{
			if (node->fRight)
				return DoFind(node->fRight, value);
			else
				return NULL;
		}
	}
	else
	{
		/* look left */
		if (value < node->fValue)
		{
			if (node->fLeft)
				return DoFind(node->fLeft, value);
			else
				return NULL;
		}
		/* look right */
		else if (value > node->fValue)
		{
			if (node->fRight)
				return DoFind(node->fRight, value);
			else
				return NULL;
		}
	}

	/* same as node->fValue */
	return node;
}

/* gather values in ascending order */
template <class TYPE>
void BinaryTreeT<TYPE>::DoAscending(BTreeNodeT<TYPE>* node, ArrayT<TYPE>& array, int& dex) const
{
	/* collect lower tree */
	if (node->fLeft) DoAscending(node->fLeft, array, dex);

	/* collect this */
	array[dex++] = node->fValue;
	
	/* collect upper tree */
	if (node->fRight) DoAscending(node->fRight, array, dex);
}

/* gaher values in descending order */
template <class TYPE>
void BinaryTreeT<TYPE>::DoDescending(BTreeNodeT<TYPE>* node, ArrayT<TYPE>& array, int& dex) const
{
	/* collect upper tree */
	if (node->fRight) DoDescending(node->fRight, array, dex);

	/* collect this */
	array[dex++] = node->fValue;

	/* collect lower tree */
	if (node->fLeft) DoDescending(node->fLeft, array, dex);
}

/* return a pointer to a new node containing the given value */
template <class TYPE>
inline BTreeNodeT<TYPE>* BinaryTreeT<TYPE>::NewNode(const TYPE& value)
{
	fSize++;
	return new BTreeNodeT<TYPE>(value);
}

/* delete the given node */
template <class TYPE>
inline void BinaryTreeT<TYPE>::FreeNode(BTreeNodeT<TYPE>* node)
{
	if (node) {
		delete node;
		fSize--;
	}
}

} // namespace Tahoe 
#endif /* _BINARY_TREE_T_H_ */
