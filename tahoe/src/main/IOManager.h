/* $Id: IOManager.h,v 1.1.1.1 2001-01-29 08:20:22 paklein Exp $ */
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
class InputBaseT;
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
	virtual void WriteOutput(double time, int ID, const dArray2DT& n_values,
		const dArray2DT& e_values);

	/* (temporarily) re-route output */
	void DivertOutput(const StringT& outfile);
	void RestoreOutput(void);

	/* accessors */
	const OutputSetT& OutputSet(int ID) const;

	// input
	int NumElementGroups(void) const;
	int NumSideSets(void) const;
	int NumNodeSets(void) const;

	void GroupNumbers(iArrayT& groupnums) const;
	void SideSetNumbers(iArrayT& sidenums) const;
	void NodeSetNumbers(iArrayT& nodenums) const;

	void ReadConnectivity(const int group, GeometryT::CodeT & geocode, iArray2DT& connects,
		iArrayT& elementmap);
	void ReadCoordinates(dArray2DT& coords, iArrayT& nodemap);
	void ReadSideSet(const int set_num, iArray2DT& sides, const int global);
	void ReadNodeSet(const int set_num, iArrayT& nodes);

	void ReadTimeSteps (dArrayT& steps);
	void ReadLabels (ArrayT<StringT>& nlabels, ArrayT<StringT>& elabels, int group_id);
	void ReadVariables (int step, int group_id, dArray2DT& nvalues, dArray2DT& evalues);

	void ReadParameters (ifstreamT& in, bool interactive, const StringT& program_name, const StringT& version);

	virtual void Interactive (void);
	virtual void ReadInputFile(ifstreamT& in);

	// use int instead of enum, so can pass derived class enum values.
	virtual void InputData (int& data, int key) const;
	virtual void InputData (iArrayT& data, int key) const;

	void SetInput(void);

	void Translate (void);

protected:

	bool ReadWord1 (ifstreamT& in, StringT& word1) const;
	virtual void Parse (ifstreamT& in, StringT& word1);
	void InteractiveIO (void);

private:

	void ReadOutputFormat (ifstreamT& in);
	void ReadInputFormat (ifstreamT& in);
	void PrintFormat (ostream& log) const;

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

	/* input formatter */
	IOBaseT::FileTypeT fInputFormat;
	InputBaseT* fInput;
	StringT     fInDatabase;

	/* echo interactive data to input file */
	ofstream fEchoInput;
	bool fEcho;
	bool fExternTahoeII;
	
private:

	/* store main out during diverions */
	OutputBaseT* fOutput_tmp;	
};

#endif
