/* $Id: iConsoleObjectT.h,v 1.5 2003-11-10 22:14:15 cjkimme Exp $ */
/* created: paklein (12/21/2000) */

#ifndef _I_CONSOLE_OBJECT_T_H_
#define _I_CONSOLE_OBJECT_T_H_

/* base class */
#include "iConsoleBaseT.h"

namespace Tahoe {

/** interface for a console object */
class iConsoleObjectT: public iConsoleBaseT
{
public:

	/** constructor */
	iConsoleObjectT(void);

	/** destructor */
	virtual ~iConsoleObjectT(void);

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

} // namespace Tahoe 
#endif /* _I_CONSOLE_OBJECT_T_H_ */
