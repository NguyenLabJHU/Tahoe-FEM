/* $Id: IOManager.cpp,v 1.9.2.1 2001-10-16 22:13:43 sawimme Exp $ */
/* created: sawimme (10/12/1999)                                          */
/* this class creates InputBaseT and OutputBaseT pointers                 */

#include "IOManager.h"

#include "fstreamT.h"
#include "ifstreamT.h"
#include "OutputSetT.h"
#include "dArrayT.h"

// output
#include "FE_ASCIIT.h"
#include "ExodusOutputT.h"
#include "EnSightOutputT.h"
#include "AbaqusOutputT.h"
#include "TecPlotOutputT.h"

IOManager::IOManager(ostream& outfile, const StringT& program_name,
	const StringT& version, const StringT& title, const StringT& input_file,
	IOBaseT::FileTypeT output_format):
	fLog(outfile),
	fOutputFormat(output_format),
	fOutput(NULL),
	fEcho (false),
	fOutputTime(0.0),
	fOutput_tmp(NULL)
{
	/* construct output formatter */
	fOutput = SetOutput(program_name, version, title, input_file, fOutputFormat);
}

IOManager::IOManager(ifstreamT& in, const IOManager& io_man):
	fLog(io_man.fLog),
	fOutputFormat(io_man.fOutputFormat),
	fOutput(NULL),
	fEcho (false),
	fOutputTime(0.0),
	fOutput_tmp(NULL)
{
	/* construct output formatter */
	fOutput = SetOutput((io_man.fOutput)->CodeName(),
				(io_man.fOutput)->Version(),
				(io_man.fOutput)->Title(), in.filename(), fOutputFormat);
}

// constructor to use in conjunction with ReadParameters
IOManager::IOManager (ostream& out) :
	fLog (out),
	fOutput(NULL),
	fOutputFormat (IOBaseT::kExodusII),
	fEcho (false)
{

}

IOManager::~IOManager(void)
{
	/* in case output is diverted */
	RestoreOutput();
	
	delete fOutput;
	fOutput = NULL;
}

void IOManager::EchoData (ostream& o) const
{
  IOBaseT temp (o);
  o << " Output format . . . . . . . . . . . . . . . . . = " << fOutputFormat  << '\n';
  temp.OutputFormats (o);
}

/*********** OUTPUT **************/


void IOManager::NextTimeSequence(int sequence_number)
{
	fOutput->NextTimeSequence(sequence_number);
}

/* set model coordinates */
void IOManager::SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_map)
{
	fOutput->SetCoordinates(coordinates, node_map);
}
	
/* register the output for an element set. returns the output ID */
int IOManager::AddElementSet(const OutputSetT& output_set)
{
	return fOutput->AddElementSet(output_set);
}

const ArrayT<OutputSetT*>& IOManager::ElementSets(void) const
{
	return fOutput->ElementSets();
}

void IOManager::AddNodeSet(const iArrayT& nodeset, int setID)
{
	fOutput->AddNodeSet(nodeset, setID);
}

void IOManager::AddSideSet(const iArray2DT& sideset, int setID, int group_ID)
{
	fOutput->AddSideSet(sideset, setID, group_ID);
}

/* output functions */
void IOManager::WriteGeometry(void)
{
	fOutput->WriteGeometry();
}

void IOManager::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT format) const
{
	if (!fOutput)
	{
		cout << "\n IOManager::WriteGeometryFile: output must be configured" << endl;
		throw eGeneralFail;		
	}

	fOutput->WriteGeometryFile(file_name, format);
}

void IOManager::WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
	fOutput->WriteOutput(fOutputTime, ID, n_values, e_values);
}

/* (temporarily) re-route output */
void IOManager::DivertOutput(const StringT& outfile)
{
	/* can only divert once */
	if (fOutput_tmp != NULL)
		cout << "\n IOManager::DivertOutput: cannot divert output to \""
		     << outfile << "\".\n"
		     <<   "     Output is already diverted to \"" <<
		     fOutput->OutputRoot() << "\"" << endl;
	else
	{
		/* store main out */
		fOutput_tmp = fOutput;
	
		/* construct temporary output formatter */
		StringT tmp(outfile);
		tmp.Append(".ext"); //OutputBaseT takes root of name passed in
		fOutput = SetOutput(fOutput_tmp->CodeName(), fOutput_tmp->Version(),
			fOutput_tmp->Title(), tmp, fOutputFormat);
		
		/* add all output sets */
		const ArrayT<OutputSetT*>& element_sets = fOutput_tmp->ElementSets();
		for (int i = 0; i < element_sets.Length(); i++)
			AddElementSet(*(element_sets[i]));
			
		/* set coordinate data */
		fOutput->SetCoordinates(fOutput_tmp->Coordinates(), fOutput_tmp->NodeMap());		
	}
}

void IOManager::RestoreOutput(void)
{
	if (fOutput_tmp != NULL)
	{
		/* delete temp output formatter */
		delete fOutput;
	
		/* restore main out */
		fOutput = fOutput_tmp;
		fOutput_tmp = NULL;
	}
}

const OutputSetT& IOManager::OutputSet(int ID) const
{
	return fOutput->OutputSet(ID);
}

/************************************************************************
* Private
************************************************************************/

/* construct and return new output formatter */
OutputBaseT* IOManager::SetOutput(const StringT& program_name,
	const StringT& version, const StringT& title, const StringT& input_file,
	IOBaseT::FileTypeT output_format)
{
	ArrayT<StringT> outstrings (4);
	outstrings[0] = input_file;
	outstrings[1] = title;
	outstrings[2] = program_name;
	outstrings[3] = version;

	const int kdigits = 4;
	OutputBaseT* output = NULL;
	switch (output_format)
	  {
	  case IOBaseT::kExodusII:
	    output = new ExodusOutputT(fLog, outstrings);
	    break;
	  case IOBaseT::kTahoe:
	  case IOBaseT::kTahoeII:
	    output = new FE_ASCIIT(fLog, true, outstrings);
	    break;
	  case IOBaseT::kEnSight:
	    output = new EnSightOutputT (fLog, outstrings, kdigits, false);
	    break;
	  case IOBaseT::kEnSightBinary:
	    output = new EnSightOutputT (fLog, outstrings, kdigits, true);
	    break;
	  case IOBaseT::kAbaqus:
	    output = new AbaqusOutputT (fLog, outstrings, false);
	    break;
	  case IOBaseT::kAbaqusBinary:
	    output = new AbaqusOutputT (fLog, outstrings, true);
	    break;
	  case IOBaseT::kTecPlot:
	    output = new TecPlotOutputT (fLog, outstrings, kdigits);
	    break;
	  default:
	    {			
	      cout << "\n IOManager::SetOutput unknown output format:"
		   << output_format << endl;
	      fLog  << "\n IOManager::SetOutput unknown output format:"
		    << output_format << endl;
	      throw eBadInputValue;
	    }
	  }	
	if (!output) throw eOutOfMemory;
	return output;
}

