// DEVELOPMENT
#include "OutPutLatticeT.h"

#include "fstreamT.h"
#include "ifstreamT.h"
#include "OutputSetT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

#include "ExceptionCodes.h"

// output
#include "FE_ASCIIT.h"
#include "ExodusOutputT.h"
#include "EnSightOutputT.h"
#include "AbaqusOutputT.h"
#include "TecPlotOutputT.h"
#include "ParaDynOutputT.h"


using namespace Tahoe;


OutPutLatticeT::OutPutLatticeT(ostream& outfile, 
			       const StringT& program_name,
			       const StringT& version, const StringT& title, 
			       const StringT& input_file,
			       IOBaseT::FileTypeT output_format, dArray2DT bounds,
			       iArrayT type):
  fLog(outfile),
  fOutputFormat(output_format),
  fOutput(NULL),
  fOutputTime(0.0)
{
  /* construct output formatter */
  fOutput = NewOutput(program_name, version, title, input_file, 
		      fOutputFormat, fLog, bounds,type);
}


OutPutLatticeT::~OutPutLatticeT(void)
{
	delete fOutput;
	fOutput = NULL;
}

/* construct and return new output formatter */
OutputBaseT* OutPutLatticeT::NewOutput(const StringT& program_name,
				       const StringT& version, 
				       const StringT& title, 
				       const StringT& input_file,
				       IOBaseT::FileTypeT output_format, 
				       ostream& log,
				       dArray2DT bounds,iArrayT type)
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
	    output = new ExodusOutputT(log, outstrings);
	    break;
	  case IOBaseT::kTahoe:
	  case IOBaseT::kTahoeII:
	  case IOBaseT::kTahoeResults:
	    output = new FE_ASCIIT(log, true, outstrings);
	    break;
	  case IOBaseT::kEnSight:
	    output = new EnSightOutputT(log, outstrings, kdigits, false);
	    break;
	  case IOBaseT::kEnSightBinary:
	    output = new EnSightOutputT(log, outstrings, kdigits, true);
	    break;
	  case IOBaseT::kAbaqus:
	    output = new AbaqusOutputT(log, outstrings, false);
	    break;
	  case IOBaseT::kAbaqusBinary:
	    output = new AbaqusOutputT(log, outstrings, true);
	    break;
	  case IOBaseT::kTecPlot:
	    output = new TecPlotOutputT(log, outstrings, kdigits);
	    break;
	  case IOBaseT::kParaDyn:
	    output = new ParaDynOutputT(log, outstrings,bounds,type);
	    break;
	  default:
	    {			
	      cout << "\n OutPutLatticeT::SetOutput unknown output format:"
		   << output_format << endl;
	      log  << "\n OutPutLatticeT::SetOutput unknown output format:"
		    << output_format << endl;
	      throw eBadInputValue;
	    }
	  }	
	if (!output) throw eOutOfMemory;
	return output;
}


void OutPutLatticeT::WriteGeometry (void)
{
  fOutput->WriteGeometry();
}

void OutPutLatticeT::WriteGeometryFile(const StringT& file_name,
	IOBaseT::FileTypeT format) const
{
	if (!fOutput)
	{
		cout << "\n OutPutLatticeT::WriteGeometryFile: output must be configured" << endl;
		throw eGeneralFail;		
	}

	fOutput->WriteGeometryFile(file_name, format);
}

void OutPutLatticeT::SetCoordinates(const dArray2DT& coordinates, const iArrayT* node_id)
{
	fOutput->SetCoordinates(coordinates, node_id); 
}



/* register the output for an element set. returns the output ID */
int OutPutLatticeT::AddElementSet(const OutputSetT& output_set)
{
	return fOutput->AddElementSet(output_set);
}

const ArrayT<OutputSetT*>& OutPutLatticeT::ElementSets(void) const
{
	return fOutput->ElementSets();
}

void OutPutLatticeT::AddNodeSet(const iArrayT& nodeset, const StringT& setID)
{
	fOutput->AddNodeSet(nodeset, setID);
}

void OutPutLatticeT::WriteOutput(int ID, const dArray2DT& n_values, const dArray2DT& e_values)
{
	fOutput->WriteOutput(fOutputTime, ID, n_values, e_values);
}

const OutputSetT& OutPutLatticeT::OutputSet(int ID) const
{
	return fOutput->OutputSet(ID);
}

