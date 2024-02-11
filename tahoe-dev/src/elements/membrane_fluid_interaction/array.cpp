#include "array.h"


namespace memFluid{


bool operator < (const array& A, const array& B){

    if(A.ptx < B.ptx){
	return true;
    }
    else if(A.ptx==B.ptx) {
	if(A.pty < B.pty){
	    return true;
	}
	else if(A.pty==B.pty) {
	    return false;
	}	
	else {
	    return false;
	}
    }
    else {
	return false;
    }

} // <

} // end memFluid
