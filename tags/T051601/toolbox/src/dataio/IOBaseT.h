/* $Id: IOBaseT.h,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: sawimme (09/28/1999)                                          */
/* Base class for InputBaseT and OutputBaseT                              */

#ifndef _IOBASE_T_H_
#define _IOBASE_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

class IOBaseT
{
public:

	enum OutputModeT {kAtFail =-2,
	                 kAtFinal =-1,
	                 kAtNever = 0,
	                   kAtInc = 1};

	enum FileTypeT {kTahoe = 0,
	              kTahoeII = 1,
	              kTecPlot = 2,
	              kEnSight = 3,
            kEnSightBinary = 4,
	             kExodusII = 5,
                   kAbaqus = 6,
             kAbaqusBinary = 7};
	
	/* constructor */
	IOBaseT(ostream& out);
	
	/* destructor */
	virtual ~IOBaseT(void);

	/* convert integer to FileTypeT */
	static FileTypeT int_to_FileTypeT(int i);
	friend istream& operator>>(istream& in, IOBaseT::FileTypeT& file_type);

protected:

	/* format the output stream */
	void SetStreamPrefs(ostream& stream) const;

	/* returns 1 if the stream is open */
	int IsOpen(ofstream& stream) const;
	
protected:
	
	ostream& fout;
};

#endif // _IOBASE_T_H_
