#ifndef GLOBFUNCS_H
#define GLOBFUNCS_H

#include "realtypes.h"
#include "parameters.h"
#include "matrix.h"
#include <list>
#include <vector>
#include <unistd.h>

namespace memFluid{

// calculate area of the node matrix
REAL polyarea(std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, 	      std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end );

int inpoly(const REAL x, const REAL y, 
std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end );	// true is 1, false is 0

REAL heaviside(const REAL val);

matrix delaunay(std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, 	          std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end );

matrix tricheck(std::vector<REAL>::const_iterator nodex_beg, 
 	        std::vector<REAL>::const_iterator nodey_beg, matrix &conn);

matrix project_point_to_line_segment(matrix &pt1, matrix &pt2, matrix &pt);

std::vector<int> sort30(std::list<REAL>::const_iterator pt_beg,
		        std::list<REAL>::const_iterator pt_end );	// return the index of the smallest 30 values

std::vector<int> convert_sort(std::vector<int>::iterator sort_beg,
		         std::vector<int>::iterator sort_end );


int sign(REAL x); // the same sign in matlab

void sparseMatrixSolver(matrix &dsol, matrix &F, int n, 
	std::vector<REAL>::iterator vals_beg, std::vector<REAL>::iterator vals_end, 
	std::vector<int>::iterator rowInds_beg, std::vector<int>::iterator rowInds_end, 
	std::vector<int>::iterator colPtr_beg, std::vector<int>::iterator colPtr_end );


size_t getTotalSystemMemory();

} // end of memFluid

#endif



