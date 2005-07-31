/* $Id: MonitorT.cpp,v 1.1.1.1 2001-01-25 20:56:25 paklein Exp $ */
/* created: paklein (09/30/1996)                                          */

#include "MonitorT.h"
#include <iomanip.h>
#include "Constants.h"

/* parameters */
const int kPercentHeadRoom = 5;

/* constructor */
MonitorT::MonitorT(void):
	fChanged(0),
	fStatusWrapper(kPercentHeadRoom, fStatus)
{

}

MonitorT::MonitorT(int dim):
	fChanged(0),
	fStatus(dim),
	fStatusWrapper(kPercentHeadRoom, fStatus)
{
	AllToON();
}

/* sets status flag for index to Marked */
void MonitorT::Mark(int index)
{
	if (fStatus[index] == kON)
	{
		fChanged       = 1;
		fStatus[index] = kMarked;		
	}
}

/* reset Marked flags */
void MonitorT::MarkedToOFF(void)
{
	fChanged = 0;
	fStatus.ChangeValue(kMarked, kOFF);
}

void MonitorT::MarkedToON(void)
{
	fChanged = 0;
	fStatus.ChangeValue(kMarked, kON);
}

/* using other status flags */
void MonitorT::SetFlag(int index, int to)
{
	if (fStatus[index] != to)
	{
		fChanged = 1;
		fStatus[index] = to;
	}
}

/* reset functions */
void MonitorT::Reset(void) { fChanged = 0; }
void MonitorT::ResetFlags(int from, int to)
{
	fChanged = 0;
	fStatus.ChangeValue(from,to);
}

/* reset all status flags to ON */
void MonitorT::AllToON(void)
{
	fStatus = kON;
}

/*
* Print indexes (real numbering) flagged as kMonitorOFF
* to the output stream with a header.
*/
void MonitorT::PrintOFF(ostream& out) const
{
	int flaggedcount = fStatus.Count(kOFF);

	out << " Number of flagged values. . . . . . . . . . . . = ";
	out << flaggedcount << "\n\n";
	
	if (flaggedcount > 0)
	{
		PrintValued(out, kOFF, 5);
		out << '\n';
	}
}

/*
* Print indexes (real numbering) with the prescribed flag value
* to the output stream.
*/
void MonitorT::PrintValued(ostream& out, int value, int wrapat) const
{
	int* p = fStatus.Pointer();
	int  count = 0;

	for (int i = 0; i < fStatus.Length(); i++)
		if (*p++ == value)
		{
			out << setw(kIntWidth) << i + 1;
			if (++count == wrapat)
			{
				out << '\n';
				count = 0;
			}	
		}
		
	/* final wrap */		
	if (count != 0)	out << '\n';
}

/* read/write all status flags */
void MonitorT::ReadStatus(istream& in)
{
	in >> fStatus;
}

void MonitorT::WriteStatus(ostream& out) const
{
	out << fStatus << '\n';
}

/* resize with the given status */
void MonitorT::Resize(int length, bool copy_in)
{
	fChanged = 1;
	fStatusWrapper.SetLength(length, copy_in);
}
