/* $Id: IOBaseT.cpp,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: sawimme (09/28/1999)                                          */

#include "IOBaseT.h"

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include "ExceptionCodes.h"

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
		default:
			cout << "\n int_to_IOFileType: could not convert: " << i << endl;
			throw eOutOfRange;
	}
	
	/* dummy */
	return IOBaseT::kTahoe;
}

istream& operator>>(istream& in, IOBaseT::FileTypeT& file_type)
{
	int i_type;
	in >> i_type;
	file_type = IOBaseT::int_to_FileTypeT(i_type);

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
