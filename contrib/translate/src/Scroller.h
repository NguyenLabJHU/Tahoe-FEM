/* $Id: Scroller.h,v 1.1 2004-11-15 09:21:19 paklein Exp $ */
#ifndef _SCROLLER_H_
#define _SCROLLER_H_

/* base class */
#include "TranslateIOManager.h"

namespace Tahoe {

/** transform ConveyorT results to scrolling view */
class Scroller: public TranslateIOManager
{
public:

	/** constructor */
	Scroller(ostream& message, istream& in, bool write);

	/** run */
	virtual void Translate(const StringT& program, const StringT& version, 
		const StringT& title);

protected:

	/** set input sources */
	void SetInput(void);

private:

	/** init output sets */
	void InitOutputSets(void);

private:

	/** \name input source */
	/*@{*/
	StringT fSource;
	IOBaseT::FileTypeT fFileType;
	/*@}*/
	
	/** output ID's */
	int fOutputID;
	
	/** \name tracking info */
	/*@{*/
	/** nodes on the cleavage plane */
	iArrayT fNodes;
	
	/** cleavage plane */
	double fCleavagePlane;
	
	/** displacement jump threshold */
	double fJumpThreshold;
	
	/** +1 for traveling right, -1 for traveling left */
	int fDirection;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _SCROLLER_H_ */
