/* $Id: ParameterUtils.h,v 1.3 2003-09-03 23:41:59 paklein Exp $ */
#ifndef _PARAMETER_UTILS_H_
#define _PARAMETER_UTILS_H_

/* direct members */
#include "ParameterInterfaceT.h"

namespace Tahoe {

/** list of integers. Defines a list of integers with default ParameterInterfaceT name 
 * "IntegerList" which contains zero or more "Integer" entries. */
class IntegerListT: public ParameterInterfaceT
{
public:
	/** \name constructors */
	/*@{*/
	IntegerListT(const StringT& name);
	IntegerListT(void);
	/*@}*/

	/** the list name */
	const StringT& ListName(void) const { return fListName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/** information about subordinate parameters */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** list name */
	StringT fListName;
};

/** list of double's. Defines a list of double's with default ParameterInterfaceT name 
 * "DoubleList" which contains zero or more "Double" entries. */
class DoubleListT: public ParameterInterfaceT
{
public:
	/** \name constructors */
	/*@{*/
	DoubleListT(const StringT& name);
	DoubleListT(void);
	/*@}*/

	/** the list name */
	const StringT& ListName(void) const { return fListName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);

	/** information about subordinate parameters */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** list name */
	StringT fListName;
};

/** named value */
template <class TYPE>
class NamedValueT: public ParameterInterfaceT
{
public:

	/** constructor */
	NamedValueT(const StringT& name);

	/** the value name */
	const StringT& ValueName(void) const { return fValueName; };

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/
	/** describe the parameters needed by the interface */
	virtual void DefineParameters(ParameterListT& list) const;

	/** accept parameter list.
	 * \param list input parameter list, which should be validated using ParameterInterfaceT::ValidateParameterList
	 *        to ensure the list conforms to the description defined by the interface. */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** value */
	TYPE fValue;

private:

	/** name */
	StringT fValueName;
};

template <class TYPE>
NamedValueT<TYPE>::NamedValueT(const StringT& name):
	ParameterInterfaceT(name)
{

}

template <class TYPE>
void NamedValueT<TYPE>::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);
	
	/* value name */
	list.AddParameter(fValueName, "name", ParameterListT::ZeroOrOnce);
	list.AddParameter(fValue, "value");
}

template <class TYPE>
void NamedValueT<TYPE>::TakeParameterList(const ParameterListT& list)
{
	/* get name */
	const ParameterT* value_name = list.Parameter("name");
	if (value_name)
		fValueName = *value_name;

	/* the value */
	list.GetParameter("value", fValue);
}

/** an integer with default ParameterInterfaceT name "Integer" */
class IntegerT: public NamedValueT<int>
{
public:

	/** \name constructors */
	/*@{*/
	IntegerT(void);
	IntegerT(const StringT& name);
	/*@{*/
};

/** a double with default ParameterInterfaceT name "Double" */
class DoubleT: public NamedValueT<double>
{
public:

	/** \name constructors */
	/*@{*/
	DoubleT(void);
	DoubleT(const StringT& name);
	/*@{*/
};

} /* namespace Tahoe */

#endif /* _PARAMETER_UTILS_H_ */
