/* $Id: IOBaseT.h,v 1.10 2002-10-20 22:36:52 paklein Exp $ */
/* created: sawimme (09/28/1999) */
#ifndef _IOBASE_T_H_
#define _IOBASE_T_H_

/* forward declarations */
#include "ios_fwd_decl.h"

namespace Tahoe {

/* forward declarations */
class InputBaseT;
class OutputBaseT;
class StringT;

/** database types and simple functions */
class IOBaseT
{
public:

	enum OutputModeT {kAtFail =-2,
	                 kAtFinal =-1,
	                 kAtNever = 0,
	                   kAtInc = 1};

	/** I/O file types */
	enum FileTypeT {
	kNone = -1,
	kTahoe = 0,
	              kTahoeII = 1,
	              kTecPlot = 2,
	              kEnSight = 3,
                kEnSightBinary = 4,
	             kExodusII = 5,
                       kAbaqus = 6,
                 kAbaqusBinary = 7,
                          kAVS = 8,
                    kAVSBinary = 9,
  	        kPatranNeutral = 10,
                 kTahoeResults = 11,
                      kParaDyn = 12 };
	
	/* constructor */
	IOBaseT(ostream& out);
	
	/* destructor */
	virtual ~IOBaseT(void);

	/** convert integer to FileTypeT */
	static FileTypeT int_to_FileTypeT(int i);
	friend istream& operator>>(istream& in, IOBaseT::FileTypeT& file_type);

	/** write list of input formats to log */
	void InputFormats (ostream &log) const;

	/** write list of output formats to log */
	void OutputFormats (ostream &log) const;
	
	/** try to guess the file format based on the file extension */
	static FileTypeT name_to_FileTypeT(const char* file_name);

	/** construct new input object. User is responsible for deleting the object.
	 * \param message stream InputBaseT will use to log messages */
	static InputBaseT* NewInput(FileTypeT format, ostream& message);

	/** construct a new output formatter */
	static OutputBaseT* NewOutput(const StringT& program_name, const StringT& version,
		const StringT& title, const StringT& input_file,
		FileTypeT output_format, ostream& log);	

protected:

	/* format the output stream */
	void SetStreamPrefs(ostream& stream) const;

	/* returns 1 if the stream is open */
	int IsOpen(ofstream& stream) const;
	
protected:
	
	ostream& fout;
};

} // namespace Tahoe 
#endif // _IOBASE_T_H_
