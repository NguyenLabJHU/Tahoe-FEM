/* $Id: garage.h,v 1.1.2.2 2003-05-04 22:13:39 paklein Exp $ */
#ifndef _GARAGE_H_
#define _GARAGE_H_

/* base class */
#include "ParameterInterfaceT.h"

/* forward declarations */
class window;

using namespace Tahoe;

class garage: public ParameterInterfaceT
{
public:

	/** constructor */
	garage(void);

	/** destructor */
	~garage(void);

	/** \name implementation of ParameterInterfaceT */
	/*@{*/
	virtual void DefineParameters(ParameterListT& list) const;
	virtual void TakeParameterList(const ParameterListT& list);

	virtual void SubNames(ArrayT<StringT>& names, ArrayT<ParameterListT::OccurrenceT>& occur, ArrayT<bool>& is_inline) const;
	virtual ParameterInterfaceT* NewSub(const StringT& list_name) const;
	/*@}*/

private:

	bool opener_;
	double length_;
	double width_;
	
	window* window_;
};

#endif /* _GARAGE_H_ */
