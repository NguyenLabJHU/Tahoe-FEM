/* $Id: ParameterUtils.cpp,v 1.4 2004-01-21 17:15:19 paklein Exp $ */
#include "ParameterUtils.h"

using namespace Tahoe;

/**********************************************************************
 * IntegerListT implementation
 **********************************************************************/

IntegerListT::IntegerListT(const StringT& name):
	NamedListT<IntegerParameterT>(name)
{

}

IntegerListT::IntegerListT(void):
	NamedListT<IntegerParameterT>("IntegerList")
{

}

/**********************************************************************
 * DoubleListT implementation
 **********************************************************************/

/* constructors */
DoubleListT::DoubleListT(const StringT& name):
	NamedListT<DoubleParameterT>(name)
{

}

DoubleListT::DoubleListT(void):
	NamedListT<DoubleParameterT>("DoubleList")
{

}

/**********************************************************************
 * StringListT implementation
 **********************************************************************/

/* constructors */
StringListT::StringListT(const StringT& name):
	NamedListT<StringParameterT>(name)
{

}

StringListT::StringListT(void):
	NamedListT<StringParameterT>("DoubleList")
{

}

/**********************************************************************
 * IntegerT implementation
 **********************************************************************/

/* constructors */
IntegerParameterT::IntegerParameterT(void):
	NamedParameterT<int>("Integer")
{
	fValue = 0;
}

IntegerParameterT::IntegerParameterT(const StringT& name):
	NamedParameterT<int>(name)
{
	fValue = 0;
}

/**********************************************************************
 * DoubleT implementation
 **********************************************************************/

/* constructors */
DoubleParameterT::DoubleParameterT(void):
	NamedParameterT<double>("Double")
{
	fValue = 0;
}

DoubleParameterT::DoubleParameterT(const StringT& name):
	NamedParameterT<double>(name)
{
	fValue = 0;
}

/**********************************************************************
 * DoubleT implementation
 **********************************************************************/

/* constructors */
StringParameterT::StringParameterT(void):
	NamedParameterT<StringT>("String")
{

}

StringParameterT::StringParameterT(const StringT& name):
	NamedParameterT<StringT>(name)
{

}
