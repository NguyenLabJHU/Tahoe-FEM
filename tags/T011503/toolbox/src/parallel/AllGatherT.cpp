/* $Id: AllGatherT.cpp,v 1.1 2002-12-05 08:25:19 paklein Exp $ */
#include "AllGatherT.h"
#include "CommunicatorT.h"

using namespace Tahoe;

/* constructor */
AllGatherT::AllGatherT(CommunicatorT& comm):
	MessageT(comm),
	fEqual(false),
	fTotal(0)
{

}

/* initialize gather data */
void AllGatherT::Initialize(int my_size)
{
	/* get counts from all */
	fCounts.Dimension(fComm.Size());
	fCounts = 0;
	
	/* collect counts from all */
	fComm.AllGather(my_size, fCounts);
	int total = fCounts.Sum();
	
	/* determine if data size is the same from all */
	fEqual = Same(fCounts);
	
	/* set displacements */
	if (fEqual)
	{
		fDisplacements.Dimension(fCounts);
		int offset = 0;
		for (int i = 0; i < fCounts.Length(); i++)
		{
			fDisplacements[i] = offset;
			offset += fCounts[i];
		}
	}
	else fDisplacements.Dimension(0);
}

void AllGatherT::AllGather(const ArrayT<double>& my_data, ArrayT<double>& gather)
{
	/* check */
	if (gather.Length() < fTotal) ExceptionT::SizeMismatch("AllGatherT::AllGather");
	
	/* equal sized or not */
	if (fEqual)
		fComm.AllGather(my_data, gather);
	else
		fComm.AllGather(fCounts, fDisplacements, my_data, gather);
}

void AllGatherT::AllGather(const ArrayT<int>& my_data, ArrayT<int>& gather)
{
	/* check */
	if (gather.Length() < fTotal) ExceptionT::SizeMismatch("AllGatherT::AllGather");

	/* equal sized or not */
	if (fEqual)
		fComm.AllGather(my_data, gather);
	else
		fComm.AllGather(fCounts, fDisplacements, my_data, gather);
}
