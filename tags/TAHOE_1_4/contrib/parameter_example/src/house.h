/* $Id: house.h,v 1.3 2003-08-14 01:22:43 paklein Exp $ */
#ifndef _HOUSE_H_
#define _HOUSE_H_

/* base class */
#include "ParameterInterfaceT.h"

/* direct members */
#include "roof.h"
#include "driveway.h"
#include "garage.h"

/* forward declarations */
class lawn;
class basement;
class room;

using namespace Tahoe;

class house: public ParameterInterfaceT
{
public:

	/** build style */
	enum style {
		undefined,
		colonial,
		ranch,
		split_level
	};

	/** constructor */
	house(void);

	/** destructor */
	~house(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	virtual void DefineSubs(SubListT& sub_list) const;
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, SubListT& sub_sub_list) const;
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** \name constructing choices */
	/*@{*/
	room* New_room(const StringT& room_name, bool throw_on_fail) const;
	basement* New_basement(const StringT& basement_name, bool throw_on_fail) const;
	/*@}*/

private:

	/** \name required, fixed simple parameters */
	/*@{*/
	/** house style as enumerated type */
	style style_;
	int zipcode_;
	/*@}*/

	/** \name required, fixed nested classes */
	/*@{*/
	roof roof_;
	driveway driveway_;
	garage garage1_;
	garage garage2_;
	/*@}*/

	/** \name variable components */
	/*@{*/
	/** optional component appears at most once */
	lawn* lawn_;

	/** optional component appears at most once which is one of a choice of basement subclasses */
	basement* basement_;

	/** optional list of component that must have at least ony entry taken from a choice of
	 * subclasses of room */
	ArrayT<room*> rooms_;
	/*@}*/
};

#endif /* _HOUSE_H_ */
