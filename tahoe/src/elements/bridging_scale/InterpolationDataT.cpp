/* $Id: InterpolationDataT.cpp,v 1.5 2005-04-22 00:29:54 paklein Exp $ */
#include "InterpolationDataT.h"
#include "iArray2DT.h"
#include "dArrayT.h"
#include "dArray2DT.h"

using namespace Tahoe;

/* transpose the given interpolation data */
void InterpolationDataT::Transpose(const InverseMapT& map, const RaggedArray2DT<int>& neighbors,
	const RaggedArray2DT<double>& neighbor_weights)
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

/* return the interpolation data in {row, column, value} (RCV) format */
void InterpolationDataT::GenerateRCV(iArrayT& r, iArrayT& c, dArrayT& v) const
{
	/* dimensions */
	int num_values = fNeighborWeights.Length();
	r.Dimension(num_values);
	c.Dimension(num_values);
	v.Dimension(num_values);

	/* global row id's */
	iArrayT rows;
	fMap.Forward(rows);

	/* collect values */
	int index = 0;
	iArrayT cols;
	dArrayT vals;
	for (int i = 0; i < fNeighborWeights.MajorDim(); i++) {
		fNeighborWeights.RowAlias(i, vals);
		fNeighbors.RowAlias(i, cols);
		for (int j = 0; j < vals.Length(); j++) {
			r[index] = rows[i];
			c[index] = cols[j];
			v[index] = vals[j];
			index++;
		}
	}
}

/* transpose the given interpolation data */
void InterpolationDataT::Transpose(const InverseMapT& map, const iArray2DT& neighbors,
	const dArray2DT& neighbor_weights)
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

/* return the interpolation data as a sparse matrix */
void InterpolationDataT::ToMatrix(iArrayT& r, iArrayT& c, dArrayT& v) const
{
	iArrayT fwd(fNeighborWeights.MajorDim());
	fMap.Forward(fwd);	

	int nv = fNeighborWeights.Length();
	r.Dimension(nv);
	c.Dimension(nv);
	v.Dimension(nv);
	
	iArrayT tmp(fNeighbors.Length(), fNeighbors.Pointer());
	c = tmp;
	v = fNeighborWeights.Pointer();

	int index = 0;	
	for (int i = 0; i < fwd.Length(); i++) {
		int dim = fNeighbors.MinorDim(i);
		for (int j = 0; j < dim; j++)
			r[index++] = fwd[i];
	}	
}
