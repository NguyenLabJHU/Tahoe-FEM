/* $Id: IOBaseT.h,v 1.16 2007-04-04 17:07:07 sawimme Exp $ */
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
class ifstreamT;

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
	kAutomatic = -2,
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
                      kParaDyn = 12,
                    kAbaqusINP = 13 };
	
	/* constructor */
	IOBaseT(ostream& out);
	
	/* destructor */
	virtual ~IOBaseT(void);

	/** convert integer to FileTypeT */
	static FileTypeT int_to_FileTypeT(int i);
	friend istream& operator>>(istream& in, IOBaseT::FileTypeT& file_type);

	/** open new stream with path defined relative to the given file. Throws exception
	 * if the file is not found */
	static void OpenRelative(ifstreamT& in, const StringT& file, const StringT& root_file);

	/** write list of input formats to log */
	static void InputFormats (ostream &log);

	/** write list of output formats to log */
	static void OutputFormats (ostream &log);
	
	/** try to guess the file format based on the file extension */
	static FileTypeT name_to_FileTypeT(const char* file_name);

	/** construct new input object. User is responsible for deleting the object.
	 * \param message stream InputBaseT will use to log messages */
	static InputBaseT* NewInput(FileTypeT format, ostream& message);

	/** construct a new output formatter */
	static OutputBaseT* NewOutput(const StringT& program_name, const StringT& version,
		const StringT& title, const StringT& output_file,
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
