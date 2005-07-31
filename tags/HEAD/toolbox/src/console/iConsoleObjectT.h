/* $Id: iConsoleObjectT.h,v 1.1.1.1 2001-01-25 20:56:28 paklein Exp $ */
/* created: paklein (12/21/2000)                                          */
/* iConsoleObjectT.h                                                      */

#ifndef _I_CONSOLE_OBJECT_T_H_
#define _I_CONSOLE_OBJECT_T_H_

/* base class */
#include "iConsoleBaseT.h"

class iConsoleObjectT: public iConsoleBaseT
{
public:

	/* constructor */
	iConsoleObjectT(void);

	/* subs control - return true if changed */
	bool iAddSub(iConsoleObjectT& sub);
	bool iDeleteSub(iConsoleObjectT& sub);

	/* accessors */
	iConsoleObjectT* iSuper(void) const;
	const ArrayT<iConsoleObjectT*>& iSubs(void) const;
	const StringT& iName(void) const;

//  protected:

	/* set name string */
	void iSetName(const StringT& name);	

private:

	/* one level up */
	iConsoleObjectT* fSuper;
	
	/* subs */
	AutoArrayT<iConsoleObjectT*> fSubs;
	
	/* data */
	StringT fName;
};

/* inlines */
inline iConsoleObjectT* iConsoleObjectT::iSuper(void) const { return fSuper; }
inline const ArrayT<iConsoleObjectT*>& iConsoleObjectT::iSubs(void) const
{
	return fSubs;
}
inline const StringT& iConsoleObjectT::iName(void) const { return fName; }

#endif /* _I_CONSOLE_OBJECT_T_H_ */
