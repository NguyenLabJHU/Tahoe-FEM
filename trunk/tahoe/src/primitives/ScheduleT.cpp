/* $Id: ScheduleT.cpp,v 1.6 2004-06-17 07:14:05 paklein Exp $ */
/* created: paklein (05/24/1996) */
#include "ScheduleT.h"
#include "ifstreamT.h"

#include <iostream.h>
#include <iomanip.h>

using namespace Tahoe;

/* constructors */
ScheduleT::ScheduleT(double value):
	fCurrentTime(0.0),
	fCurrentValue(value)
{

}

ScheduleT::ScheduleT(int numpts):
	fTime(numpts),
	fValue(numpts),
	fCurrentTime(0.0),
	fCurrentValue(0.0)
{

}

ScheduleT::ScheduleT(const dArrayT& times, const dArrayT& values):
	fTime(times),
	fValue(values),
	fCurrentTime(0.0),
	fCurrentValue(0.0)
{
	/* check data */
	CheckSequential();
}

/* I/O operators */
void ScheduleT::Read(ifstreamT& in)
{
	for (int i = 0; i < fTime.Length(); i++)
	{
		in >> fTime[i];
		in >> fValue[i];
	}

	/* check data */
	CheckSequential();
}

void ScheduleT::Write(ostream& out) const
{
	int d_width = out.precision() + kDoubleExtra;

	out << setw(d_width) << "time";
	out << setw(d_width) << "factor" << '\n';

	for (int i = 0; i < fTime.Length(); i++)
	{
		out << setw(d_width) << fTime[i];
		out << setw(d_width) << fValue[i] << '\n';
	}
}

/* set the load factor based on the time given */
void ScheduleT::SetTime(double time)
{
	fCurrentTime = time;
	int num_pts = fTime.Length();
	if (num_pts > 1)
	{
		if ( time <= fTime[0] )			/* first abscissa */
			fCurrentValue = fValue[0];
		else if ( time >= fTime[num_pts-1] )	/* last abscissa  */
			fCurrentValue = fValue[num_pts-1];
		else
			for (int i = 0; i < num_pts; i++)
				if (fTime[i] >= time)
				{		
					fCurrentValue = fValue[i-1] + (time - fTime[i-1])*
										(fValue[i] - fValue[i-1])/
										(fTime[i] - fTime[i-1]);
					i = num_pts;		/* exit */
				}
	}
	else if (num_pts == 1)
		fCurrentValue = fValue[0];
}

/* check that times are sequential */
void ScheduleT::CheckSequential(void) const
{
	if (fTime.Length() == 1) return;
	for (int i = 1; i < fTime.Length(); i++)
		if (fTime[i] < fTime[i-1]) throw ExceptionT::kGeneralFail;
}
