/* $Id: MapIteratorT.h,v 1.1 2008-02-11 13:57:06 paklein Exp $ */
#ifndef _MAP_ITERATOR_T_H_
#define _MAP_ITERATOR_T_H_

/* direct members */
#include "MapT.h"
#include "MapNodeT.h"
#include "AutoArrayT.h"

namespace Tahoe {

/** iterator for the map class */
template <class key_TYPE, class value_TYPE>
class MapIteratorT: protected BinaryTreeT<MapNodeT<key_TYPE, value_TYPE> >
{
public:

	/** constructor */
	MapIteratorT(MapT<key_TYPE, value_TYPE>& map);

	/** restart the iteration */
	void Top(void);

	/** return the next {key, value} in the map. Returns 1 if the next pair exists, otherwise 0  */
	int Next(key_TYPE* key, value_TYPE* value);

	/** return the current stack depth */
	int Depth(void) { return fStack.Length(); }
	int MaxDepth(void) { return fMaxDepth; }

private:

	/** the map */
	MapT<key_TYPE, value_TYPE>& fMap;
	
	/** traversal stack */
	AutoArrayT<const BTreeNodeT<MapNodeT<key_TYPE, value_TYPE> >* > fStack;
	int fMaxDepth;
	
};

/* constructor */
template <class key_TYPE, class value_TYPE>
MapIteratorT<key_TYPE, value_TYPE>::MapIteratorT(MapT<key_TYPE, value_TYPE>& map): 
	fMap(map),
	fMaxDepth(0) {
	Top();
}

/* restart the iteration */
template <class key_TYPE, class value_TYPE>
void MapIteratorT<key_TYPE, value_TYPE>::Top(void) { 
	
	/* clear out history */
	fStack.Dimension(0);
	fMaxDepth = 0;
	
	/* restart stack */
	BinaryTreeT<MapNodeT<key_TYPE, value_TYPE> >& btree = fMap.BTree();
	if (btree.Root()) {
		fStack.Push(btree.Root());
	}
};

/* next in tree */
template <class key_TYPE, class value_TYPE>
int MapIteratorT<key_TYPE, value_TYPE>::Next(key_TYPE* key, value_TYPE* value) {
	if (fStack.Length() == 0) {
		return 0;
	} else {
	
		/* get the leading node/remove from stack */
		const BTreeNodeT<MapNodeT<key_TYPE, value_TYPE> >* next = fStack[0];
		fStack.DeleteAt(0); 
		
		/* resolve values */
		const MapNodeT<key_TYPE, value_TYPE>& key_value = next->Value();
		*key   = key_value.Key();
		*value = key_value.Value();
		
		/* push left node, if not NULL */
		if (next->Left()) fStack.Push(next->Left());

		/* push right node, if not NULL */
		if (next->Right()) fStack.Push(next->Right());
	
		/* record max depth */
		fMaxDepth = (fStack.Length() > fMaxDepth) ? fStack.Length() : fMaxDepth;
	
		return 1;
	}
}

} /* namespace Tahoe */

#endif /* _MAP_ITERATOR_T_H_ */
