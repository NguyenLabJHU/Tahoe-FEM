/* $Id: InverseMapT.cpp,v 1.5 2004-04-23 20:22:31 paklein Exp $ */
#include "InverseMapT.h"
#include "iArrayT.h"

using namespace Tahoe;

/* construct the inverse map */
void InverseMapT::SetMap(const nArrayT<int>& forward)
{
	if (forward.Length() == 0)
	{
		fShift = 0;
		Dimension(0);
	}
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
		int* inv_map = Pointer();
		int dim = forward.Length();
		for (int i = 0; i < dim; i++)
			inv_map[forward[i] - fShift] = i;
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
