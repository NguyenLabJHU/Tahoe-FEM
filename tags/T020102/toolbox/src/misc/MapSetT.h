/* $Id: MapSetT.h,v 1.1 2002-01-09 18:25:25 paklein Exp $ */
/* created: paklein (03/26/2000)                                          */

#ifndef _MAP_SET_T_H_
#define _MAP_SET_T_H_

/* direct members */
#include "iArrayT.h"

class MapSetT
{
public:

	/* constructor */
	MapSetT(void) { }

	/* dimensioning */
	void Allocate(int n_sets, int e_sets) {
		fNodalMaps.Allocate(n_sets);
		fElementMaps.Allocate(e_sets);
	}
		
	/* number of maps */
	int NumNodeMaps(void) const { return fNodalMaps.Length(); }
	int NumElementMaps(void) const { return fElementMaps.Length(); }
		
	/* the maps */
	iArrayT& NodeMap(int set) { return fNodalMaps[set]; }
	iArrayT& ElementMap(int set) { return fElementMaps[set]; }
	const iArrayT& NodeMap(int set) const { return fNodalMaps[set]; }
	const iArrayT& ElementMap(int set) const { return fElementMaps[set]; }

	 private:
	
	ArrayT<iArrayT> fNodalMaps;
	ArrayT<iArrayT> fElementMaps;
};

#endif /* _MAP_SET_T_H_ */
