#include "cell.h"
#include <iostream>
#include <cstdlib>

namespace dem{

void cell::setNodes(int m, int n, int i, int j){
	if(m<=0 || n<=0 || i<=0 || j<=0){
		std::cout << "Error when setNodes for cell: m, n, i, j should be positive!" << std::endl;
		exit(-1);
	}
	if(m==n || m==i || m==j || n==i || n==j || i==j){
		std::cout << "Error when setNodes for cell: m, n, i, j should be different!" << std::endl;
		exit(-1);	
	}
	node_m = m;
	node_n = n;
	node_i = i;
	node_j = j;
}

// non member functions
std::ostream& operator << (std::ostream& out, const cell& A){
	out << A.node_m << "\t" << A.node_n << "\t" << A.node_i << "\t" << A.node_j;
	return out;
}

} // namespace dem
