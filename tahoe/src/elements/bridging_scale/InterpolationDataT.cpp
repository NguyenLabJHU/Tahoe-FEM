/* $Id: InterpolationDataT.cpp,v 1.2 2004-03-04 08:54:20 paklein Exp $ */
#include "InterpolationDataT.h"

using namespace Tahoe;

/* transpose the given InterpolationDataT */
void InterpolationDataT::Transpose(const InverseMapT& map, RaggedArray2DT<int>& neighbors,
	RaggedArray2DT<double>& neighbor_weights)
{
	/* determine unique neighbor numbers */
	iArrayT neighbors_all(neighbors.Length(), neighbors.Pointer());
	iArrayT neighbors_used;
	neighbors_used.Union(neighbors_all);
	
	/* map to data for each neighbor */
	fMap.SetMap(neighbors_used);

	/* counts */
	iArrayT neigh;
	iArrayT count(neighbors_used.Length());
	count = 0;
	for (int i = 0; i < neighbors.MajorDim(); i++) {
		neighbors.RowAlias(i, neigh);
		for (int j = 0; j < neigh.Length(); j++) {
			int col_map = fMap.Map(neigh[j]);
			count[col_map]++;
		}
	}

	/* dimension arrays */
	fNeighbors.Configure(count);
	fNeighborWeights.Configure(count);

	/* get the forward map */
	iArrayT forward;
	map.Forward(forward);

	/* transpose the tables */
	count = 0;
	dArrayT weight;
	for (int i = 0; i < neighbors.MajorDim(); i++) {
		int row = forward[i];
		neighbors.RowAlias(i, neigh);
		neighbor_weights.RowAlias(i, weight);
		for (int j = 0; j < neigh.Length(); j++) {
			int col_map = fMap.Map(neigh[j]);
			int& dex = count[col_map];
			fNeighbors(col_map, dex) = row;
			fNeighborWeights(col_map, dex) = weight[j];
			dex++;
		}
	}	
}
