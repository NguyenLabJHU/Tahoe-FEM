/* $Id: enum2int.h,v 1.3 2003-05-06 23:37:53 cjkimme Exp $ */
#ifndef _ENUM_2_INT_H_
#define _ENUM_2_INT_H_

namespace Tahoe {

/** little wrapper to let enum's be treated as int's. Works only if
 * enum's are actually treated as integers. */
template <class enum_TYPE>
class enum2int
{
public:

	/** constructor */
	enum2int(enum_TYPE& e);

	/** enum as int (lvalue) */
	operator int&() { return e_int; };

private:

	/** the enum instance */
	enum_TYPE& e_;
	
	/** enum as int */
	int& e_int;
};

template <class enum_TYPE>
enum2int<enum_TYPE>::enum2int(enum_TYPE& e):
	e_(e),
	e_int((int&) e_)
{

}

} /* namespace Tahoe */

#endif /* _ENUM_PROXY_H_ */


