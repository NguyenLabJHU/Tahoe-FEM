/* $Id: IOManager.h,v 1.5.2.1 2001-10-16 22:13:44 sawimme Exp $ */
/* created: sawimme (10/12/1999)                                          */

#ifndef _IOMANAGER_H_
#define _IOMANAGER_H_

/* language support */
#include <iostream.h>
#include <fstream.h>

/* direct members */
#include "StringT.h"
#include "iArrayT.h"
#include "iAutoArrayT.h"
#include "IOBaseT.h"
#include "GeometryT.h"

/* forward declarations */
class ifstreamT;
class OutputBaseT;
class iArray2DT;
class dArray2DT;
class dArrayT;
class OutputSetT;

class IOManager
{
public:

	/* constructor */
	IOManager(ostream& outfile, const StringT& program_name, const StringT& version,
		const StringT& title, const StringT& input_file, IOBaseT::FileTypeT output_format);
	IOManager(ifstreamT& in, const IOManager& io_man);
	IOManager (ostream& out);

	/* destructor */
	virtual ~IOManager(void);

	/* echo format file type */
	void EchoData (ostream& o) const;

	void NextTimeSequence(int sequence_number);

	/* set model coordinates */
	void SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_map);
	
	/* register the output for an element set. returns the output ID */
	int AddElementSet(const OutputSetT& output_set);
	const ArrayT<OutputSetT*>& ElementSets(void) const;
	void AddNodeSet(const iArrayT& nodeset, int setID);
	void AddSideSet(const iArray2DT& sideset, int setID, int group_ID); // ID from AddElementSet

	/* output functions */
	void WriteGeometry(void);
	void WriteGeometryFile(const StringT& file_name, IOBaseT::FileTypeT format) const;
	
	void SetOutputTime(double time);
	virtual void WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values);

	/* (temporarily) re-route output */
	void DivertOutput(const StringT& outfile);
	void RestoreOutput(void);

	/* accessors */
	const OutputSetT& OutputSet(int ID) const;

private:

	/* return new output formatter */
	OutputBaseT* SetOutput(const StringT& program_name, const StringT& version,
		const StringT& title, const StringT& input_file,
		IOBaseT::FileTypeT output_format);	

protected:

	ostream& fLog;
	StringT fTitle;

	/* output formatter */
	IOBaseT::FileTypeT fOutputFormat;
	OutputBaseT* fOutput;

	/* echo interactive data to input file */
	ofstream fEchoInput;
	bool fEcho;
	
private:

	/* run time */
	double fOutputTime;

	/* store main out during diverions */
	OutputBaseT* fOutput_tmp;	
};

/* inlines */
inline void IOManager::SetOutputTime(double time) { fOutputTime = time; }

#endif
