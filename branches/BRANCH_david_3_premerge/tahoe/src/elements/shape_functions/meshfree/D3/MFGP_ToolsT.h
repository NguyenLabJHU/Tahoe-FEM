/* $Id: MFGP_ToolsT.h,kyonten*/
/* created: kyonten */
#ifndef _MFGP_TOOLS_T_H_
#define _MFGP_TOOLS_T_H_

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"

namespace Tahoe {

/* forward declarations */
class LocalArrayT;

/* class contains functions needed for "meshfree_grad_plast" element */
class MFGP_ToolsT
{
	public:
	/* constructor */
	MFGP_ToolsT(void);
	
	/* destructor */
	~MFGP_ToolsT(void);

	/** 1D/2D/3D unsymmetric matrix (arising from the third derivatives
	 *  of the shape function) stored as index vector: resolve components
	 *  Note: out of 27 (8 in 2D) components of the matrix, only 9 (4 in 2D) are needed.
	 **/
	static void ExpandIndex3(int nsd, int dex, int& dex_1, int& dex_2, int& dex_3);
	
    static void ExpandIndex2(int nsd, int dex_1, int dex_2, int dex_3, int& dex_12, 
                             int& dex_23, int& dex_31);
    /** compute the jacobian of the nodal values.
	 * uses externally provided shape function derivatives.
	 * \param nodal values at the nodes: [nnd] x [nu]
	 * \param DDNa shape function second derivatives: [nstr] x [nnd]
	 * \param jacobian resulting jacobian: [nu] x [nstr]
	 **/
	//id JacobianD2(const LocalArrayT& nodal, const dArray2DT& DDNa, dMatrixT& jacobian) const;
	
	/** compute the jacobian of the nodal values.
	 * uses externally provided shape function derivatives.
	 * \param nodal values at the nodes: [nnd] x [nu]
	 * \param DDDNa shape function third derivatives: [nsd*nsd] x [nnd] 
	 * \param jacobian resulting jacobian: [nu] x [nsd*nsd]  
	 **/
	//id JacobianD3(const LocalArrayT& nodal, const dArray2DT& DDDNa, dMatrixT& jacobian) const;
	
};

} /* namespace Tahoe */
#endif /* _MFGP_TOOLS_T_H_ */
