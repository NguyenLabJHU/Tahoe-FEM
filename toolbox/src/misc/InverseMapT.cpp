/* $Id: InverseMapT.cpp,v 1.6 2005-02-24 17:34:03 paklein Exp $ */
#include "InverseMapT.h"
#include "iArrayT.h"

using namespace Tahoe;

/* construct the inverse map */
void InverseMapT::SetMap(const nArrayT<int>& forward)
{
	if (forward.Length() == 0)
		Free();
	else
	{
		/* range */
		int max;
		forward.MinMax(fShift, max);
		int range = max - fShift + 1;
	
		/* dimension */
		Dimension(range);
		AutoArrayT<int>::operator=(-1);

		/* make map */
		fEntrees = 0;
		int* inv_map = Pointer();
		int dim = forward.Length();
		for (int i = 0; i < dim; i++)
		{
			/* foward map must be unique */
			int& entry = inv_map[forward[i] - fShift];
			if (entry != -1)
				ExceptionT::GeneralFail("InverseMapT::SetMap", 
					"forward map contains repeated entry %d at %d", forward[i], i+1);
			else
				entry = i;
		}
		fEntrees = forward.Length();
	}
}

/* recover the forward map */
void InverseMapT::Forward(ArrayT<int>& forward) const
{
	iArrayT tmp;
	tmp.Alias(*this);
	int dim = tmp.Length() - tmp.Count(-1);
	forward.Dimension(dim);

	int dex = 0;
	for (int i = 0; i < tmp.Length(); i++)
		if (tmp[i] != -1)
			forward[dex++] = i + fShift;
}
