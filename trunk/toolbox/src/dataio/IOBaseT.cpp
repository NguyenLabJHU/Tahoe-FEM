/* $Id: IOBaseT.cpp,v 1.8 2002-07-23 11:33:42 sawimme Exp $ */
/* created: sawimme (09/28/1999) */

#include "IOBaseT.h"

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "ExceptionCodes.h"
#include "StringT.h"


using namespace Tahoe;

IOBaseT::IOBaseT(ostream& out): fout(out) { }
IOBaseT::~IOBaseT(void) { }

/* convert integer to FileTypeT */
IOBaseT::FileTypeT IOBaseT::int_to_FileTypeT(int i)
{
	switch (i)
	{
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
		default:
			cout << "\n int_to_IOFileType: could not convert: " << i << endl;
			throw eOutOfRange;
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

void IOBaseT::InputFormats (ostream& log) const
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

void IOBaseT::OutputFormats (ostream& log) const
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
}

/* try to guess the file format based on the file extension */
IOBaseT::FileTypeT IOBaseT::name_to_FileTypeT(const char* file_name)
{
	StringT ext;
	ext.Suffix(file_name);
	
	if (ext == ".exo" || ext == ".e")
		return kExodusII;
	else if (ext == ".geom")
		return kTahoeII;
	else if (ext == ".case")
		return kEnSight;
	else if (ext == ".run" || ext == ".geo")
		return kTahoeResults;
	else {
		cout << "\n IOBaseT::name_to_FileTypeT: could not guess file type from name: " 
		    << file_name << endl;
		throw eGeneralFail;
	}

	/* dummy */
	return kTahoe;
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
