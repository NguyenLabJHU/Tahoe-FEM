// file: IOBaseT.cpp

// created      : SAW (09/28/1999)
// last modified: PAK (11/07/1999)

#include "IOBaseT.h"

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "ExceptionCodes.h"

IOBaseT::IOBaseT(ostream& out): fout(out) { }
IOBaseT::~IOBaseT(void) { }

/* convert integer to IOFileType */
IOBaseT::IOFileType int_to_IOFileType(int i)
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
		default:
			cout << "\n int_to_IOFileType: could not convert: " << i << endl;
			throw eOutOfRange;
	}
	
	/* dummy */
	return IOBaseT::kTahoe;
}

istream& operator>>(istream& in, IOBaseT::IOFileType& file_type)
{
	int i_type;
	in >> i_type;
	file_type = int_to_IOFileType(i_type);

	return in;
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
