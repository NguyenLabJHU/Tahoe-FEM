/* $Id: MF_dMatrixT.h,kyonten*/
/* created: kyonten */
#ifndef _MF_D_MATRIX_T_H_
#define _MF_D_MATRIX_T_H_

namespace Tahoe {

/** interface for a 1D/2D/3D unsymmetric matrix stored as 
 * an index vector */
/* Note: out of 27 components of the matrix, only 9 are needed.  kyonten */

class MF_dMatrixT
{
public:

	/* constructor */
	MF_dMatrixT(void);
	
	/* destructor */
	~MF_dMatrixT(void);

	static void ExpandIndex3(int nsd, int dex, int& dex_1, int& dex_2, int& dex_3);
    static void ExpandIndex2(int dex_1, int dex_2, int& dex);
    
};


/* inlines */
		
inline void MF_dMatrixT::ExpandIndex3(int nsd, int dex, int& dex_1, int& dex_2, int& dex_3)
{
#if __option(extended_errorcheck)
	/* consistency check */
	const char caller[] = "MF_dMatrixT::ExpandIndex3";
	if (dex >= nsd*nsd) ExceptionT::OutOfRange(caller, "bad index %d", dex);
#endif
	
	int  map_1D[3] = {0,0,0};
	int  map_2D[12] = {0,0,0,0,1,0,1,0,0,1,1,0};
	int map_3D[27] = {0,0,0,1,1,0,2,2,0,0,0,1,1,1,1,2,2,1,
	                  0,0,2,1,1,2,2,2,2};
	int* map_list[4] = {NULL, map_1D, map_2D, map_3D};
	int* map = map_list[nsd];
	int* p = map + 3*dex; 
	dex_1 = p[0];
	dex_2 = p[1];
	dex_3 = p[2];
}


inline void MF_dMatrixT::ExpandIndex2(int dex_1, int dex_2, int& dex)
{	
	if (dex_1==0 && dex_2==0)
	{
		dex = 0;
	}
	if ((dex_1==0 || dex_1==1) && (dex_2==1 || dex_2==0))
	{
		dex = 1;
	}
	if ((dex_1==0 || dex_1==2) && (dex_2==2 || dex_2==0))
	{
		dex = 3;
	}
	if (dex_1==1 && dex_2==1)
	{
		dex = 4;
	}
	if ((dex_1==1 || dex_1==2) && (dex_2==2 || dex_2==1))
	{
		dex = 5;
	}
	if (dex_1==2 && dex_2==2)
	{
		dex = 6;
	}
}


} /* namespace Tahoe */
#endif /* _MF_D_MATRIX_T_H_ */
