/* $Id: iConsoleObjectT.h,v 1.2 2001-11-28 22:05:44 paklein Exp $ */
/* created: paklein (12/21/2000) */

#ifndef _I_CONSOLE_OBJECT_T_H_
#define _I_CONSOLE_OBJECT_T_H_

/* base class */
#include "iConsoleBaseT.h"

/** interface for a console object */
class iConsoleObjectT: public iConsoleBaseT
{
public:

	/** constructor */
	iConsoleObjectT(void);

	/** add a sub console.
	 * \return true if successfully added, false otherwise */
	bool iAddSub(iConsoleObjectT& sub);

	/** remove a sub console.
	 * \return true if successfully removed, false otherwise */
	bool iDeleteSub(iConsoleObjectT& sub);

	/** return a pointer to the console super */
	iConsoleObjectT* iSuper(void) const;
	const ArrayT<iConsoleObjectT*>& iSubs(void) const;
	const StringT& iName(void) const;

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
