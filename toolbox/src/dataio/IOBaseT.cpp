/* $Id: IOBaseT.cpp,v 1.17 2005-04-28 20:58:53 paklein Exp $ */
/* created: sawimme (09/28/1999) */
#include "IOBaseT.h"

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "ExceptionT.h"
#include "StringT.h"

/* input formats */
#include "TahoeInputT.h"
#include "ExodusInputT.h"
#include "PatranInputT.h"
#include "EnSightInputT.h"
#include "AbaqusInputT.h"
#include "TextInputT.h"

/* output formats */
#include "TextOutputT.h"
#include "ExodusOutputT.h"
#include "EnSightOutputT.h"
#include "AbaqusOutputT.h"
#include "TecPlotOutputT.h"
#include "ParaDynOutputT.h"

using namespace Tahoe;

IOBaseT::IOBaseT(ostream& out): fout(out) { }
IOBaseT::~IOBaseT(void) { }

/* convert integer to FileTypeT */
IOBaseT::FileTypeT IOBaseT::int_to_FileTypeT(int i)
{
	switch (i)
	{
		case -2:
			return IOBaseT::kAutomatic;
		case -1:
			return IOBaseT::kNone;
		case 0:
			return IOBaseT::kTahoe;
		case 1:
			return IOBaseT::kTahoeII;
		case 2:			
			return IOBaseT::kTecPlot;
		case 3:
			return IOBaseT::kEnSight;
		case 4:
			return IOBaseT::kEnSightBinary;
		case 5:
			return IOBaseT::kExodusII;
		case 6:
			return IOBaseT::kAbaqus;
		case 7:
			return IOBaseT::kAbaqusBinary;
   	        case 8:
	                return IOBaseT::kAVS;
	        case 9:
	                return IOBaseT::kAVSBinary;
	       case 10:
	                return IOBaseT::kPatranNeutral;
	       case 11:
	                return IOBaseT::kTahoeResults;
	       case 12:
	                return IOBaseT::kParaDyn;
		default:
			ExceptionT::OutOfRange("IOBaseT::int_to_FileTypeT", 
				"could not convert %d", i);
	}
	
	/* dummy */
	return IOBaseT::kTahoe;
}

namespace Tahoe {

istream& operator>>(istream& in, IOBaseT::FileTypeT& file_type)
{
	int i_type;
	in >> i_type;
	file_type = IOBaseT::int_to_FileTypeT(i_type);

	return in;
}

}

/* open new stream with path defined relative to the given file */
void IOBaseT::OpenRelative(ifstreamT& in, const StringT& file, const StringT& root_file)
{
	/* build file path */
	StringT filename(file);
	filename.ToNativePathName();
	StringT path;
	path.FilePath(root_file);
	filename.Prepend(path);

	/* open stream */
	in.open(filename);
	if (!in.is_open())
		ExceptionT::DatabaseFail("IOBaseT::OpenRelative", 
			"could not open file \"%s\"", filename.Pointer());
}

void IOBaseT::InputFormats (ostream& log)
{
  log << "    eq. " << setw (2) << IOBaseT::kTahoe         << ". Tahoe\n";
  log << "    eq. " << setw (2) << IOBaseT::kTahoeII       << ". Tahoe II\n";
//log << "    eq. " << setw (2) << IOBaseT::kTecPlot       << ". TecPlot 7.5\n";
  log << "    eq. " << setw (2) << IOBaseT::kEnSight       << ". Ensight 6 Gold ASCII\n";
  log << "    eq. " << setw (2) << IOBaseT::kEnSightBinary << ". Ensight 6 Gold Binary\n";
  log << "    eq. " << setw (2) << IOBaseT::kExodusII      << ". Exodus II\n";
  log << "    eq. " << setw (2) << IOBaseT::kAbaqus        << ". ABAQUS ASCII (.fin)\n";
  log << "    eq. " << setw (2) << IOBaseT::kAbaqusBinary  << ". ABAQUS Binary (.fil)\n";
//log << "    eq. " << setw (2) << IOBaseT::kAVS           << ". AVS UCD ASCII\n";
//log << "    eq. " << setw (2) << IOBaseT::kAVSBinary     << ". AVS UCD Binary\n";
  log << "    eq. " << setw (2) << IOBaseT::kPatranNeutral << ". PATRAN Neutral\n";
  log << "    eq. " << setw (2) << IOBaseT::kTahoeResults  << ". Tahoe Results (.geo/.run)\n";
}

void IOBaseT::OutputFormats (ostream& log)
{
//log << "    eq. " << setw (2) << IOBaseT::kTahoe         << ". Tahoe\n";
  log << "    eq. " << setw (2) << IOBaseT::kTahoeII       << ". Tahoe II\n";
  log << "    eq. " << setw (2) << IOBaseT::kTecPlot       << ". TecPlot 7.5\n";
  log << "    eq. " << setw (2) << IOBaseT::kEnSight       << ". Ensight 6 Gold ASCII\n";
  log << "    eq. " << setw (2) << IOBaseT::kEnSightBinary << ". Ensight 6 Gold Binary\n";
  log << "    eq. " << setw (2) << IOBaseT::kExodusII      << ". Exodus II\n";
  log << "    eq. " << setw (2) << IOBaseT::kAbaqus        << ". ABAQUS ASCII (.fin)\n";
  log << "    eq. " << setw (2) << IOBaseT::kAbaqusBinary  << ". ABAQUS Binary (.fil)\n";
  log << "    eq. " << setw (2) << IOBaseT::kAVS           << ". AVS UCD ASCII\n";
//log << "    eq. " << setw (2) << IOBaseT::kAVSBinary     << ". AVS UCD Binary\n";
  log << "    eq. " << setw (2) << IOBaseT::kPatranNeutral << ". PATRAN Neutral\n";
  log << "    eq. " << setw (2) << IOBaseT::kParaDyn       << ". PARADYN\n";
}

/* try to guess the file format based on the file extension */
IOBaseT::FileTypeT IOBaseT::name_to_FileTypeT(const char* file_name)
{
	StringT ext;
	ext.Suffix(file_name);
	
	if (ext == ".exo" || ext == ".e" || ext == ".g" || ext == ".gen")
		return kExodusII;
	else if (ext == ".geom")
		return kTahoeII;
	else if (ext == ".case")
		return kEnSight;
	else if (ext == ".run" || ext == ".geo")
		return kTahoeResults;
	else if (ext == ".atoms")
		return kParaDyn;
	else
		ExceptionT::GeneralFail("IOBaseT::name_to_FileTypeT",
			"could not guess file type from \"%s\"", file_name);

	/* dummy */
	return kTahoe;
}

/* construct new input object */
InputBaseT* IOBaseT::NewInput(FileTypeT format, ostream& message)
{
	const char caller[] = "IOBaseT::NewInput";
	InputBaseT* input = NULL;
	try {
	switch (format)
    {
		case kTahoe:
      		/* do nothing, arrays will be registered via ElementBaseT and NodeManager */
			input = NULL;
			break;

		case kAutomatic: /* default type */
		case kTahoeII:
			input = new TahoeInputT(message);
			break;

		case kEnSight:
			input = new EnSightInputT(message, false);
      		break;

		case kEnSightBinary:
			input = new EnSightInputT(message, true);
			break;

#ifdef __ACCESS__
		case kExodusII:
			input = new ExodusInputT(message);
			break;
#endif

		case kPatranNeutral:
			input = new PatranInputT(message);
			break;

		case kAbaqus:
		case kAbaqusBinary:
			input = new AbaqusInputT(message);
			break;

		case kTahoeResults:
			input = new TextInputT(message);
			break;

		default:
			ExceptionT::GeneralFail(caller, "unsupported module format %d", format);
    }
    } /* end try */
    
    catch(ExceptionT::CodeT exception) {
    	ExceptionT::Throw(exception, caller);
    }    
    return input;
}

/* construct and return new output formatter */
OutputBaseT* IOBaseT::NewOutput(const StringT& program_name,
	const StringT& version, const StringT& title, const StringT& output_file,
	FileTypeT output_format, ostream& log)
{
	const char caller[] = "IOBaseT::NewOutput";
	ArrayT<StringT> outstrings (4);
	outstrings[0] = output_file;
	outstrings[1] = title;
	outstrings[2] = program_name;
	outstrings[3] = version;

	const int kdigits = 4;
	OutputBaseT* output = NULL;
	try {
	switch (output_format)
	  {
	  case IOBaseT::kExodusII:
	    output = new ExodusOutputT(log, outstrings);
	    break;
	  case IOBaseT::kAutomatic:
	  case IOBaseT::kTahoe:
	  case IOBaseT::kTahoeII:
	  case IOBaseT::kTahoeResults:
	    output = new TextOutputT(log, true, outstrings);
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
	    output = new ParaDynOutputT(log, outstrings);
	    break;
	  default:
	  	ExceptionT::GeneralFail(caller, "unknown output format %d", output_format);
	}
	} /* end try */  

    catch(ExceptionT::CodeT exception) {
    	ExceptionT::Throw(exception, caller);
    }    
	return output;
}

/*************************************************************************
* Protected
*************************************************************************/

/* format the output stream */
void IOBaseT::SetStreamPrefs(ostream& stream) const
{
	stream.precision(12);
	stream.setf(ios::showpoint);
	stream.setf(ios::right, ios::adjustfield);
	stream.setf(ios::scientific, ios::floatfield);
}

/* returns 1 if the stream is open */
int IOBaseT::IsOpen(ofstream& stream) const
{
#ifdef __MWERKS__
	return stream.is_open();
#else
// is_open is only defined for filebuf not ostream or istream,
// and isn't defined as const
ofstream* non_const_ifstr = (ofstream*) &stream;
filebuf* fbuf = non_const_ifstr->rdbuf();
return fbuf->is_open();
#endif
}
