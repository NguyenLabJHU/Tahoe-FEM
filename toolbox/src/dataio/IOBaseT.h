/* $Id: IOBaseT.h,v 1.2 2001-06-14 12:55:10 sawimme Exp $ */
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
                 kAbaqusBinary = 7,
                          kAVS = 8,
                    kAVSBinary = 9,
  	        kPatranNeutral = 10 };
	
	/* constructor */
	IOBaseT(ostream& out);
	
	/* destructor */
	virtual ~IOBaseT(void);

	/* convert integer to FileTypeT */
	static FileTypeT int_to_FileTypeT(int i);
	friend istream& operator>>(istream& in, IOBaseT::FileTypeT& file_type);

	void PrintFormat (ostream &log) const;

protected:

	/* format the output stream */
	void SetStreamPrefs(ostream& stream) const;

	/* returns 1 if the stream is open */
	int IsOpen(ofstream& stream) const;
	
protected:
	
	ostream& fout;
};

#endif // _IOBASE_T_H_
