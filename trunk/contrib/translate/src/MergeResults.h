/* $Id: MergeResults.h,v 1.2 2002-07-02 21:22:59 cjkimme Exp $ */
#ifndef _TRANSLATE_MERGE_H_
#define _TRANSLATE_MERGE_H_

/* base class */
#include "TranslateIOManager.h"

namespace Tahoe {

/** merge results data from multiple results files */
class MergeResults: public TranslateIOManager
{
public:

	/** constructor */
	MergeResults (ostream& message, istream& in, bool write);

	/** destructor */
	virtual ~MergeResults(void);

	/** run */
	virtual void Translate (const StringT& program, const StringT& version, 
		const StringT& title);

protected:

	/** set input sources */
	void SetInput(void);

private:

	/** generate combined coordinates list
	 * \param coords combined, expanded coordinate list. The rows in this
	 *        row correspond to the (global) id of the node.
	 * \param nodes_used list nodes in the expanded list that were used */
	void CombinedCoordinates(dArray2DT& coords, iArrayT& nodes_used);
		
private:

	/** array of input sources */
	ArrayT<ModelManagerT*> fInputs;

	/** node maps from the input sources */
	ArrayT<iArrayT> fNodeMaps;
};

} // namespace Tahoe

#endif /* _TRANSLATE_MERGE_H_ */
