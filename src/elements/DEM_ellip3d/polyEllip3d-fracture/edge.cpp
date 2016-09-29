#include "edge.h"
#include <iostream>
#include <cstdlib>

namespace dem{

void edge::setNodes(int m, int n){
	if(m<=0 || n<=0){
		std::cout << "Error when setNodes for edge: m, n should be positive!" << std::endl;
		exit(-1);
	}
	if(m==n){
		std::cout << "Error when setNodes for edge: m, n should be different!" << std::endl;
		exit(-1);	
	}
	if(m>n){	// swap m, n to make sure m < n
		int temp = m;
		m = n;
		n = temp;
	}

	node_m = m;
	node_n = n;
}

// non member functions
/*inline */bool operator < (const edge &A, const edge &B){
	if (A.node_m==0 || A.node_n==0 || B.node_m==0 || B.node_n==0){
		std::cout << "Error when compare two objects of edge: nodes of edges should be positive!" << std::endl;
		exit(-1);
	}	
	if(A.node_m < B.node_m){
		return true;
	}
	else if(A.node_m==B.node_m) {
		if(A.node_n < B.node_n){
			return true;
		}
		else if(A.node_n==B.node_n) {
			return false;
		}	
		else {
			return false;
		}
	}
	else {
		return false;
	}
}

std::ostream& operator << (std::ostream& out, const edge& A){
	out << A.node_m << "\t" << A.node_n;
	return out;
}

} // namespace dem

/*
// used for debug
using namespace dem;
int main() {
	edge e1;
	std::cout << "e1: " << std::endl;
	std::cout << e1 << std::endl;
	e1.setNodes(11,19);
	std::cout << "e1: " << std::endl;
	std::cout << e1 << std::endl;
	e1.setNodes(13,11);
	std::cout << "e1: " << std::endl;
	std::cout << e1 << std::endl;
	edge e2;
	e2.setNodes(13,11);
	if(e1 < e2)
		std::cout << "e1<e2" << std::endl;
	else
		std::cout << "e1>e2" << std::endl;
}
*/





