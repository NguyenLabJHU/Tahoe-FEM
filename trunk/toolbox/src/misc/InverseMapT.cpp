/* $Id: InverseMapT.cpp,v 1.1 2002-11-25 07:08:30 paklein Exp $ */
#include "InverseMapT.h"
#include "nArrayT.h"

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
