#include "globfuncs.h"
#include "sparse.h"
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace memFluid{


// calculate area of the node matrix
// refer to http://mathworld.wolfram.com/PolygonArea.html
REAL polyarea(std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, 	      std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end){

    std::vector<REAL> nodex(nodex_beg, nodex_end);
    std::vector<REAL> nodey(nodey_beg, nodey_end);

    REAL area = 0;
    std::vector<REAL>::size_type n = nodex.size();
    for (std::vector<REAL>::size_type i=0; i<n-1; i++){

        std::vector<REAL>::size_type j = (i + 1);
//        area += (nodey[i] + nodey[j])*(nodex[i] - nodex[j]);
	area += nodex[i]*nodey[j]-nodex[j]*nodey[i];
    }
    area += nodex[n-1]*nodey[0]-nodex[0]*nodey[n-1];	// XnY1-X1Yn

    return area*0.5;

// Polygon is CCW if area > 0.0f
// Polygon is CW if area < 0.0f
} // polyarea()



matrix tricheck(std::vector<REAL>::const_iterator nodex_beg, 
 	        std::vector<REAL>::const_iterator nodey_beg, matrix &conn){

    int count = 0;

    matrix tri = conn;

    for(int e=1; e<conn.num_row+1; ++e){

	matrix sctr;
	sctr.appendRow(conn.getRow(e)); // conn is nx3, sctr is 1x3
	matrix node_sctr(3,2);	// node_sctr is 3x2
	std::vector<REAL>::const_iterator itx1, ity1;
	itx1 = nodex_beg+sctr(1,1)-1;
	ity1 = nodey_beg+sctr(1,1)-1;	// sctr is numbering starting from 1
	std::vector<REAL>::const_iterator itx2, ity2;
	itx2 = nodex_beg+sctr(1,2)-1;
	ity2 = nodey_beg+sctr(1,2)-1;
	std::vector<REAL>::const_iterator itx3, ity3;
	itx3 = nodex_beg+sctr(1,3)-1;
	ity3 = nodey_beg+sctr(1,3)-1;
	node_sctr(1,1) = (*itx1); node_sctr(1,2) = (*ity1);
	node_sctr(2,1) = (*itx2); node_sctr(2,2) = (*ity2);
	node_sctr(3,1) = (*itx3); node_sctr(3,2) = (*ity3);

	matrix N(3,1);
	REAL xi = 1.0/3.0;
	REAL eta = 1.0/3.0;
	N(1,1) = 1-xi-eta; N(2,1) = xi; N(3,1) = eta;
	matrix dNdxi(3,2);
	dNdxi(1,1) = -1; dNdxi(1,2) = -1;
	dNdxi(2,1) = 1;  dNdxi(2,2) = 0;
	dNdxi(3,1) = 0;  dNdxi(3,2) = 1;

	matrix temp = node_sctr.getTrans()*dNdxi;
	REAL detJ = det(temp);	// 2x2

	if(detJ<0){
	    tri(e,1) = sctr(1,3); tri(e,2) = sctr(1,2); tri(e,3) = sctr(1,1);
	    count = count+1;
	}
	else if(detJ==0){
	    std::cout << "ZERO JACOBI IN ELEMENT " << e << " CANNOT FIX" << std::endl;
	}
    }

    return tri;

} // tricheck()


int inpoly(const REAL x, const REAL y, 
std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end ){
// true is 1, false is 0

    REAL TOL = 1.0e-12;
    
    std::vector<REAL> nodex(nodex_beg, nodex_end);
    std::vector<REAL> nodey(nodey_beg, nodey_end);

    int nnode = nodex.size();

    std::vector<int> edge_x;	// the x nodes of the edges, starts from 1
    std::vector<int> edge_y;

    for(int i=1; i<nnode; ++i ){
	edge_x.push_back(i);
	edge_y.push_back(i+1);
    }
    edge_x.push_back(nnode);
    edge_y.push_back(1);


//// PRE-PROCESSING
//  int  n = 1; // only one point needed to be checked
    int nc = edge_x.size();

    REAL dxymin = 100;
    REAL tol = TOL*dxymin;

    // MAIN LOOP
    bool cn = false;

    bool on = cn;
    for(std::vector<int>::size_type k=1; k<nc+1; ++k){

	// nodes in current edge
	int n1 = edge_x[k-1];
	int n2 = edge_y[k-1];

	// Endpoints
	REAL y1 = nodey[n1-1];
	REAL y2 = nodey[n2-1];

	REAL x1, x2;
	if(y1<y2){
	    x1 = nodex[n1-1];
	    x2 = nodex[n2-1];
	}
	else{
	    REAL yt = y1;
	    y1 = y2;
	    y2 = yt;
	    x1 = nodex[n2-1];
	    x2 = nodex[n1-1];
	}
	REAL xmin, xmax;
	if(x1>x2){
	    xmin = x2;
	    xmax = x1;
	}
	else{
	    xmin = x1;
	    xmax = x2;
	}

	// Binary search to find first point with y<=y1 for current edge
	// no needs, since we only have one point

	// check point
	if(y<=y2){
	    if(x>=xmin){
		if(x<=xmax){
		    // check if we are "on" the edge
		    if( fabs( (y2-y)*(x1-x)-(y1-y)*(x2-x) ) < tol )
			on = true;
		
		    // do the actual intersection test
		    if( y<y2 && (y2-y1)*(x-x1)<(y-y1)*(x2-x1) )
			cn = !cn;
		} // xmax
	    } // xmin
	    else if (y<y2){
		cn = !cn;
	    }
	} // y2

    } // end for

    cn = cn||on;

    int in;
    if(cn==true)
	in = 1;	// 1 means true
    else
	in = 0;

    return in;

} // inpoly()


REAL heaviside(const REAL X){

    REAL val;
    if(X<0)
	val = 0;
    else if(X==0)
	val = 0.5;
    else
	val = 1;

    return val;

} // heaviside()


matrix project_point_to_line_segment(matrix&A, matrix&B, matrix&p){
// returns q, the closest point to p on the line segment from A to B

    matrix q;
    // vector from A to B
    matrix AB = B-A;
    // squared distance from A to B
    REAL AB_squared = AB(1,1)*AB(1,1)+AB(1,2)*AB(1,2);
    if(AB_squared==0){
  	// A and B are the same point
	q = A;
    }
    else{
	// vector from A to p
        matrix Ap = p-A;
    	// from http://stackoverflow.com/questions/849211/
    	// Consider the line extending the segment, parameterized as A + t (B - A)
    	// We find projection of point p onto the line.
    	// It falls where t = [(p-A) . (B-A)] / |B-A|^2

        REAL t = (Ap(1,1)*AB(1,1)+Ap(1,2)*AB(1,2))/AB_squared;
	if(t < 0.0){
	    // "Before" A on the line, just return A
	    q = A;
	}
	else if(t > 1.0){
	    // "After" B on the line, just return B
	    q = B;	
	}
	else{
	    // projection lines "inbetween" A and B on the line
	    q = A + t*AB;
	}

    } // end else

    return q;

} // project_point_to_line_segement()



std::vector<int> sort30(std::list<REAL>::const_iterator pt_beg,
		        std::list<REAL>::const_iterator pt_end ){
// return the index of the smallest 30 values
// index starts from 1

    std::vector<std::pair<REAL, int> > vp;
    int index=1;	// index
    for(std::list<REAL>::const_iterator it=pt_beg; it!=pt_end; ++it, ++index){
	vp.push_back(std::make_pair((*it), index));
    }

    // Sorting will put lower values ahead of larger ones,
    // resolving ties using the original index
    std::sort(vp.begin(), vp.end());

    std::vector<int> index_vec;
    int num = 0;
    for(std::vector<std::pair<REAL, int> >::const_iterator it=vp.begin(); it<vp.end(); ++it, ++num){
	index_vec.push_back((*it).second);
	if(num>30)
	    break;
    }

    return index_vec;

} // sort30()


std::vector<int> convert_sort(std::vector<int>::iterator sort_beg,
		         std::vector<int>::iterator sort_end ){

    std::vector<int> index_vec(sort_beg, sort_end);
    // convert 2&3, 4&5, 6&7, 8&9
    for(std::vector<int>::iterator it=index_vec.begin()+1; it<index_vec.end(); ++it){
	if((it+1)!=index_vec.end()){
	    int temp = (*it);	// temp is 2nd
	    (*it) = *(it+1);	// it is 2nd
	    it++;		// it is 3rd
	    (*it) = temp;
	}
    }

    return index_vec;

} // convert_sort()


matrix delaunay(std::vector<REAL>::const_iterator nodex_beg, std::vector<REAL>::const_iterator nodex_end, 	          std::vector<REAL>::const_iterator nodey_beg, std::vector<REAL>::const_iterator nodey_end ){
    
    std::vector<REAL>::size_type num = nodex_end-nodex_beg;

    //create input for Qhull
    std::ofstream ofs("input_for_Qhull");
    if(!ofs){
	std::cout << "Stream error when create input for Qhull!" << std::endl;
	exit (-1);
    }
    ofs.setf(std::ios::scientific, std::ios::floatfield);
    ofs.precision(OPREC);
    ofs << std::setw(OWID) << "2" << std::endl
	<< std::setw(OWID) << num << std::endl;

    std::vector<REAL>::const_iterator ity = nodey_beg;
    for(std::vector<REAL>::const_iterator itx=nodex_beg; itx<nodex_end; ++itx, ++ity){

    	ofs << std::setw(16) << (*itx)
	    << std::setw(16) << (*ity) << std::endl;
    }

    ofs.close();

    // call qhull
    //std::cout << "Checking if processor is available..." << std::endl;
    //if(system(NULL)) std::cout << "Ok!" << std::endl;
    //else exit(EXIT_FAILURE);
    system("./qdelaunay Qt i < input_for_Qhull TO tess_info");	// call the external command qdelaunay

    // read tessellation
    std::ifstream ifs("tess_info");
    if(!ifs) {
	std::cout << "stream error!" << std::endl; exit(-1);
    }
    int m, n, i;
    int totalNum;	// total number of cells, ie triangles

    ifs>>totalNum;

    matrix tri(totalNum, 3);
    for(int it=1; it!=totalNum+1; it++){
	ifs>>m>>n>>i;	// read nodes in each line
	m = m+1;
	n = n+1;
	i = i+1;	//the ID from Qhull is starting from 0

   	tri(it,1) = m;
	tri(it,2) = n;
	tri(it,3) = i;
    }

    return tri;

} // delaunay()


int sign(REAL x){
 
    if(x>0)
	return 1;
    if(x==0)
	return 0;
    if(x<0)
	return -1;

} // sign()


void sparseMatrixSolver(matrix &dsol, matrix &F, int n, 
	std::vector<REAL>::iterator vals_beg, std::vector<REAL>::iterator vals_end, 
	std::vector<int>::iterator rowInds_beg, std::vector<int>::iterator rowInds_end, 
	std::vector<int>::iterator colPtr_beg, std::vector<int>::iterator colPtr_end ) {

    int num_vals = vals_end - vals_beg;
    int num_row  = rowInds_end - rowInds_beg;
    int num_col = colPtr_end - colPtr_beg;
	
    double* vals = new double[num_vals];
    int* rowInds = new int[num_row];
    int* colPtr = new int[num_col];

    int i=0;
    for(std::vector<REAL>::iterator it=vals_beg; it<vals_end; ++it, i)
	vals[i] = (*it);

    i=0;
    for(std::vector<int>::iterator it=rowInds_beg; it<rowInds_end; ++it, i)
	rowInds[i] = (*it);

    i=0;
    for(std::vector<int>::iterator it=colPtr_beg; it<colPtr_end; ++it, i)
	colPtr[i] = (*it);	


    SparseMatrix* s = new SparseMatrix(n, n, vals, rowInds, colPtr);

    SparseMatrix* decomposed = s->cholesky_decompose_lower_triangular_returning_lower_triangular();
    delete s;

    // get b, namely F
//    double* b = new double[n];
//    for(int ii=0; ii<n; ++ii)
//	b[ii] = F(ii+1, 1);	// F is column vector

    std::vector<REAL> b;
    for(int ii=1; ii<n+1; ++ii)
	b.push_back( F(ii,1) );

    double* x = decomposed->solve_eqn(b.begin(), b.end());
    delete decomposed;

    // get dsol
    for(int ii=0; ii<n; ++ii)
	dsol(ii+1, 1) = x[ii];

    delete[] x;

} // sparseMatrixSolver()


size_t getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
} // getTotalSystemMemory()



} // end of memFluid
