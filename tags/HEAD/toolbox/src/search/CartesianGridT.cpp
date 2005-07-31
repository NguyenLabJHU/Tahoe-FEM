/* $Id: CartesianGridT.cpp,v 1.1.1.1 2001-01-25 20:56:26 paklein Exp $ */
/* created: paklein (11/10/2000)                                          */

#include "CartesianGridT.h"
#include <math.h>
#include "GraphT.h"

static int to_integer(double a) { return int((2.0*a + 1.0)/2.0); };

/* constructor */
CartesianGridT::CartesianGridT(void)
{

}

/* (re-) set dimensions */
void CartesianGridT::SetDimensions(const iArrayT& dimensions,
	const ArrayT<BoundaryConditionT>& bc)
{
	/* store */
	fDimensions = dimensions;
	fBC = bc;
	
	/* check */
	if (dimensions.Length() != fBC.Length())
	{
		cout << "\n CartesianGridT::SetDimensions: dimensionality " << dimensions.Length()
		     << " does not\n" <<   "     match the number of boundary conditions "
		     << fBC.Length() << endl;
		throw eSizeMismatch;
	}	

	/* indexed distances */
	int n_dim = fDimensions.Length();
	fShift.Allocate(n_dim);
	int shift = 1;
	for (int i = n_dim - 1; i > -1; i--)
	{
		fShift[i] = shift;	
		shift *= fDimensions[i];
	}	
	fWidth = fShift;
	for (int j = 0; j < n_dim; j++)
		fWidth[j] *= fDimensions[j] - 1;
}	

/* partition the grid (num_parts >= num_cells) */
void CartesianGridT::PartitionGrid(int num_parts, const iArrayT& cell_weight)
{
	/* checks */
	int num_cells = fDimensions.Product();
	if (num_parts > num_cells)
	{
		cout << "\n CartesianGridT::PartitionGrid: number of partitions " << num_parts
		     << " must be <= number\n"
		     <<   "     of cells in the grid " << num_cells << endl;
		throw eGeneralFail;
	}
	if (cell_weight.Length() != num_cells)
	{
		cout << "\n CartesianGridT::Partition: number of weights "
		     << cell_weight.Length() << " does not\n"
		     <<   "     match the number of grid cells " << num_cells << endl;
		throw eSizeMismatch;
	}

	/* allocate space */
	fCellMap.Allocate(num_cells);
	int num_neighbors = to_integer(pow(2, fDimensions.Length()));
	fNeighborList.Allocate(num_cells, num_neighbors + 1);
	
	/* set neighbors data */
	SetNeighborLists();
	
	iArray2DT* to_graph= &fNeighborList;
	
#if 0
	/* reduce to pair list */
	int no_neigh = fNeighborList.Count(-1);
	iArray2DT pairs(num_cells*num_neighbors - no_neigh, 2);
	int pair = 0;
	for (int i = 0; i < num_cells; i++)
	{
		int* list = fNeighborList(i);
		for (int j = 0; j < num_neighbors; j++)
			if (list[j] != -1)
			{
				int* edge = pairs(pair++);
				edge[0] = i;
				edge[1] = list[j];
			}
	}
	if (pair != pairs.MajorDim()) throw eGeneralFail;
	to_graph= &pairs;
#endif
	
	/* quick exit */
	if (num_cells == num_parts)
		fCellMap.SetValueToPosition();
	else /* load balance */
	{
		/* initialize */
		fCellMap = -1;
		
		/* initialize graph object */
		GraphT graph;
		graph.AddGroup(*to_graph);
	
		/* make graph */
		graph.MakeGraph();
		
		/* partition and balance */
		iArrayT config(1);
		config = num_parts;
		graph.Partition(config, cell_weight, fCellMap, true);
	}
}

/**********************************************************************
* Private
**********************************************************************/

/* recursive function to set neighbors */
void CartesianGridT::SetNeighborLists(void)
{
	/* index array */
	int num_dims = fDimensions.Length();
	iArrayT index(num_dims);
	index = 0;
	
	int cell_num = 0;
	while (cell_num < fNeighborList.MajorDim())
	{
		/* assign neighbors on inner-most index*/
		for (int& i = index.Last(); i < fDimensions.Last(); i++)
		{
			int* neighbors = fNeighborList(cell_num);
			for (int j = 0; j < num_dims; j++)
			{
				if (index[j] == 0)
				{
					if (fBC[j] == kFree)
						neighbors[0] = -1;
					else if (fBC[j] == kPeriodic)
						neighbors[0] = cell_num + Width(j);
					else
						throw eGeneralFail;

					neighbors[1] = cell_num + Shift(j);
				}
				else if (index[j] == fDimensions[j] - 1)
				{
					neighbors[0] = cell_num - Shift(j);

					if (fBC[j] == kFree)
						neighbors[1] = -1;
					else if (fBC[j] == kPeriodic)
						neighbors[1] = cell_num - Width(j);
					else
						throw eGeneralFail;
				}
				else
				{
					neighbors[0] = cell_num - Shift(j);
					neighbors[1] = cell_num + Shift(j);
				}
				
				/* next */
				neighbors += 2;
			}
			
			/* add self */
			*neighbors = cell_num;
			
			/* next */
			cell_num++;
		}
	
		/* increment indices */
		for (int k = num_dims - 1; k > -1; k--)
			if (index[k] == fDimensions[k])
			{
				index[k] = 0;
				if (k > 0) index[k-1]++;
			}		
	}
}
