/* $Id: MapNodeT.h,v 1.3 2003-03-09 20:36:41 paklein Exp $ */
#ifndef _MAP_NODE_T_H_
#define _MAP_NODE_T_H_

/* direct members */
#include "ExceptionT.h"

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
	
	/** default constructor. Requires key_TYPE to have a default constructor  */
	MapNodeT(void);
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

	/** \name accessors */
	/*@{*/
	const key_TYPE& Key(void) const { return fKey; };
	
	/** value must be set otherwise an exception will be thrown */
	const value_TYPE& Value(void) const;
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

/* default constructor */
template <class key_TYPE, class value_TYPE>
inline MapNodeT<key_TYPE, value_TYPE>::MapNodeT(void):
	fValue(NULL)
{

}

/* destructor */
template <class key_TYPE, class value_TYPE>
inline MapNodeT<key_TYPE, value_TYPE>::~MapNodeT(void)
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
inline bool MapNodeT<key_TYPE, value_TYPE>::operator<(const MapNodeT& rhs) const
{
	return this->fKey < rhs.fKey;
}

/* true if (the key of this) == (the key of rhs) */
template <class key_TYPE, class value_TYPE>
inline bool MapNodeT<key_TYPE, value_TYPE>::operator==(const MapNodeT& rhs) const
{
	return this->fKey == rhs.fKey;
}

/* value must be set otherwise an exception will be thrown */
template <class key_TYPE, class value_TYPE>
inline const value_TYPE& MapNodeT<key_TYPE, value_TYPE>::Value(void) const
{
	if (!fValue) ExceptionT::GeneralFail("MapNodeT::Value", "value not set");
	return *fValue;
}

} /* namespace Tahoe */

#endif /* _MAP_NODE_T_H_ */
