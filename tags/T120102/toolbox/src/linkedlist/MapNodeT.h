/* $Id: MapNodeT.h,v 1.1 2002-11-16 20:44:59 paklein Exp $ */
#ifndef _MAP_NODE_T_H_
#define _MAP_NODE_T_H_

namespace Tahoe {

/** forward declaration */
template <class key_TYPE, class value_TYPE> class MapT;

/** container MapT
 * Implementation of operator needed by MapT such that comparisons
 * of MapNodeT's depend only on the key values. */
template <class key_TYPE, class value_TYPE>
class MapNodeT
{
	friend class MapT<key_TYPE, value_TYPE>;

public:

	/** \name constructors */
	/*@{*/
	MapNodeT(const key_TYPE &key, const value_TYPE &value);

	/** copy constructor */
	MapNodeT(const MapNodeT& source);

	/** node with key only */
	MapNodeT(const key_TYPE &key);
	/*@}*/

	/** destructor */
	~MapNodeT(void);

	/** assignment operator */
	MapNodeT& operator=(const MapNodeT& rhs);

	/** \name operators for comparing key values */
	/*@{*/
	/** true if (the key of this) > (the key of rhs) */
	bool operator>(const MapNodeT& rhs) const;

	/** true if (the key of this) < (the key of rhs) */
	bool operator<(const MapNodeT& rhs) const;

	/** true if (the key of this) == (the key of rhs) */
	bool operator==(const MapNodeT& rhs) const;
	/*@}*/
	 	 	
private:

	/** \name key-value pair */
	/*@{*/
	key_TYPE    fKey;
	value_TYPE* fValue;
	/*@}*/
};

/*************************************************************************
 * Implementation
 *************************************************************************/

/* constructor */
template <class key_TYPE, class value_TYPE>
inline MapNodeT<key_TYPE, value_TYPE>::MapNodeT(const key_TYPE &key, const value_TYPE &value):
	fKey(key),
	fValue(NULL)
{
	fValue = new value_TYPE(value);
}

/** copy constructor */
template <class key_TYPE, class value_TYPE>
inline MapNodeT<key_TYPE, value_TYPE>::MapNodeT(const MapNodeT& source):
	fKey(source.fKey),
	fValue(NULL)
{
	fValue = new value_TYPE(*source.fValue);
}

/* node with key only */
template <class key_TYPE, class value_TYPE>
inline MapNodeT<key_TYPE, value_TYPE>::MapNodeT(const key_TYPE &key):
	fKey(key),
	fValue(NULL)
{

}

/* destructor */
template <class key_TYPE, class value_TYPE>
MapNodeT<key_TYPE, value_TYPE>::~MapNodeT(void)
{
	delete fValue;
}

/* assignment operator */
template <class key_TYPE, class value_TYPE>
inline MapNodeT<key_TYPE, value_TYPE>& MapNodeT<key_TYPE, value_TYPE>::operator=(const MapNodeT& rhs)
{
	fKey = rhs.fKey;
	if (fValue) delete fValue;
	fValue = new value_TYPE(*rhs.fValue);
	return *this;
}

/* true if (the key of this) > (the key of rhs) */
template <class key_TYPE, class value_TYPE>
inline bool MapNodeT<key_TYPE, value_TYPE>::operator>(const MapNodeT& rhs) const
{
	return this->fKey > rhs.fKey;
}

/* true if (the key of this) < (the key of rhs) */
template <class key_TYPE, class value_TYPE>
bool MapNodeT<key_TYPE, value_TYPE>::operator<(const MapNodeT& rhs) const
{
	return this->fKey < rhs.fKey;
}

/* true if (the key of this) == (the key of rhs) */
template <class key_TYPE, class value_TYPE>
bool MapNodeT<key_TYPE, value_TYPE>::operator==(const MapNodeT& rhs) const
{
	return this->fKey == rhs.fKey;
}

} /* namespace Tahoe */

#endif /* _MAP_NODE_T_H_ */
