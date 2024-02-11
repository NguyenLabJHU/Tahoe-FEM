#ifndef EDGE_H
#define EDGE_H

#include <iostream>

namespace dem{

class edge{
      private:
	int node_m, node_n;
      public:
	edge():node_m(0), node_n(0) {}	// constructor
	void setNodes(int, int);	// set values for node_m and node_n, make sure node_m<node_n
	int getm() const {return node_m;}
	int getn() const {return node_n;}
//	void clear();	// clear the values of node_m and node_n to zeros

	// non member functions
	friend /*inline */bool operator < (const edge&, const edge&);	// compare two edges based on m&n
	friend std::ostream& operator << (std::ostream&, const edge&);	// overload << for output

}; 

} // namespace dem



#endif
