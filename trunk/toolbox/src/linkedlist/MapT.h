/* $Id: MapT.h,v 1.1 2002-11-16 20:44:59 paklein Exp $ */
#ifndef _MAP_T_H_
#define _MAP_T_H_

/* base class */
#include "BinaryTreeT.h"

/* direct members */
#include "MapNodeT.h"

namespace Tahoe {

/** map class.
 * create a list in which keys of an arbitrary type can refer to
 * values of a second type. The key_TYPE must satisfy the requirements
 * needed to within a BinaryTreeT. In addition, the key type requires
 * operator<< */
template <class key_TYPE, class value_TYPE>
class MapT: protected BinaryTreeT<MapNodeT<key_TYPE, value_TYPE> >
{
public:

	/** constructor */
	MapT(void) {};
	
	/** return the size of the map */
	int Size(void) { return BinaryTreeT<MapNodeT<key_TYPE, value_TYPE> >::Size(); };
	
	/** insert value in the map. Returns true if the value was added, returns false if
	 * the key value already exists. */
	bool Insert(const key_TYPE& key, const value_TYPE& value);

	/** read/write access to values */
	value_TYPE& operator[](const key_TYPE& key);
};

/* insert value in the map */
template <class key_TYPE, class value_TYPE>
bool MapT<key_TYPE, value_TYPE>::Insert(const key_TYPE& key, const value_TYPE& value)
{
	MapNodeT<key_TYPE, value_TYPE> node(key, value);
	return InsertUnique(node);
}

/* read/write access to values */
template <class key_TYPE, class value_TYPE>
value_TYPE& MapT<key_TYPE, value_TYPE>::operator[](const key_TYPE& key)
{
	/* return the tree node with the given value or NULL if not present */
	MapNodeT<key_TYPE, value_TYPE> find_node(key);
	BTreeNodeT<MapNodeT<key_TYPE, value_TYPE> >* tree_node = Find(find_node);
	
	/* no match */
	if (!tree_node) {
		cout << "\n MapT<key_TYPE, value_TYPE>::operator[]: key node found: " << key << endl;
		throw ExceptionT::kOutOfRange;
	}
	
	/* return value */
	const MapNodeT<key_TYPE, value_TYPE>& node = tree_node->Value();
	return *(node.fValue);
}

} /* namespace Tahoe */

#endif /* _MAP_T_H_ */
