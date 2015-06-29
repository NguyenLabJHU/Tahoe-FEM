#ifndef CELL_H
#define CELL_H

#include <iostream>
#include "realtypes.h"
#include "matrix.h"

namespace dem{

class cell{
      private:
	int node_m, node_n, node_i, node_j;	// nodes of this cell, ie pyramid
	REAL cellVolume_init;			// the initial cell volume, used for lagrangian granular strain
	matrix bigB_init;			// the binitial bigB matrix, used for lagrangian granular strain
      public:
	cell():node_m(0), node_n(0), node_i(0), node_j(0) {} // constructor
	void setNodes(int, int, int, int);	// set values for m,n,i,j, no need to make m<n<i<j
	void setInitialCellVolume(REAL v) {cellVolume_init = v;}
	void setInitialBigB(matrix b) {bigB_init = b;}
	int getm() const {return node_m;}
	int getn() const {return node_n;}
	int geti() const {return node_i;}
	int getj() const {return node_j;}
	REAL getInitialCellVolume() const {return cellVolume_init;}
	matrix getInitialBigB() const {return bigB_init;}
//	void clear();	// clear the values of m,n,i,j to be zeros
	
	// non member functions
	friend std::ostream& operator << (std::ostream& out, const cell&);	// overload << for output

};



} // namespace dem

#endif
