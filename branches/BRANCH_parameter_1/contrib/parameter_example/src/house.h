/* $Id: house.h,v 1.1.2.3 2003-05-03 09:08:27 paklein Exp $ */
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
	virtual void SetParameters(const ParameterListT& list);

	virtual void SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur,
		ArrayT<bool>& is_inline) const;
	virtual void DefineInlineSub(const StringT& sub, ParameterListT::ListOrderT& order, 
		ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur, ArrayT<bool>& is_inline) const;
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	/** \name constructing choices */
	/*@{*/
	room* New_room(const StringT& room_name, bool throw_on_fail) const;
	basement* New_basement(const StringT& basement_name, bool throw_on_fail) const;
	/*@}*/

private:

	/** \name fixed components */
	/*@{*/
	/** house style as enumerated type */
	style style_;
	int zipcode_;
	roof roof_;
	driveway driveway_;
	garage garage1_;
	garage garage2_;
	/*@}*/

	/** \name variable components */
	/*@{*/
	ArrayT<room*> rooms_;
	lawn* lawn_;
	basement* basement_;
	/*@}*/
};

#endif /* _HOUSE_H_ */
