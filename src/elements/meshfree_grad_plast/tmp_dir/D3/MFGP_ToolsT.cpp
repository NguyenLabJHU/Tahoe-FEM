/* $Id: MFGP_ToolsT.cpp,kyonten*/
/* created: kyonten */
#include "MFGP_ToolsT.h"
//#include "toolboxConstants.h"
#include "dArray2DT.h"
#include "LocalArrayT.h"
#include "dSymMatrixT.h"
//#include "ParentDomainT.h"
using namespace Tahoe;


namespace Tahoe {
DEFINE_TEMPLATE_STATIC const bool ArrayT<MFGP_ToolsT*>::fByteCopy = true;
DEFINE_TEMPLATE_STATIC const bool ArrayT<const MFGP_ToolsT*>::fByteCopy = true; 
DEFINE_TEMPLATE_STATIC const bool ArrayT<MFGP_ToolsT>::fByteCopy = false; 
}

void MFGP_ToolsT::ExpandIndex3(int nsd, int dex, int& dex_1, int& dex_2, int& dex_3)
{
#if __option(extended_errorcheck)
	/* consistency check */
	const char caller[] = "MFGP_ToolsT::ExpandIndex3";
	if (dex >= nsd*nsd) ExceptionT::OutOfRange(caller, "bad index %d", dex);
#endif
	
	int  map_1D[3] = {0,0,0};
	int  map_2D[12] = {0,0,0,1,1,0,0,0,1,1,1,1};
	int map_3D[27] = {0,0,0,1,1,0,2,2,0,
	                  0,0,1,1,1,1,2,2,1,
	                  0,0,2,1,1,2,2,2,2};
	int* map_list[4] = {NULL, map_1D, map_2D, map_3D};
	int* map = map_list[nsd];
	int* p = map + 3*dex; 
	dex_1 = p[0];
	dex_2 = p[1];
	dex_3 = p[2];
}


void MFGP_ToolsT::ExpandIndex2(int nsd, int dex_1, int dex_2, int dex_3, 
            int& dex_12, int& dex_23, int& dex_31)
{
#if __option(extended_errorcheck)
	/* consistency check */
	const char caller[] = "MFGP_ToolsT::ExpandIndex2";
#endif	
	switch (nsd)
	{
		case 1:
		{
		 if (dex_1==0 && dex_2==0 && dex_3==0)
		 {
		 	dex_12=0; dex_23=0; dex_31=0;
		 }
		 break;	
		}
		
		case 2:
		{
			if (dex_1==0 && dex_2==0 && dex_3==0)
		 	{
		 		dex_12=0; dex_23=0; dex_31=0;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==0)
		 	{
		 		dex_12=1; dex_23=2; dex_31=2;
		 	}
		 	if (dex_1==0 && dex_2==0 && dex_3==1)
		 	{
		 		dex_12=0; dex_23=2; dex_31=2;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==1)
		 	{
		 		dex_12=1; dex_23=1; dex_31=1;
		 	}
		 	break;
		}
		
		case 3:
		{
			if (dex_1==0 && dex_2==0 && dex_3==0)
		 	{
		 		dex_12=0; dex_23=0; dex_31=0;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==0)
		 	{
		 		dex_12=1; dex_23=5; dex_31=5;
		 	}
		 	if (dex_1==2 && dex_2==2 && dex_3==0)
		 	{
		 		dex_12=2; dex_23=4; dex_31=4;
		 	}
		 	if (dex_1==0 && dex_2==0 && dex_3==1)
		 	{
		 		dex_12=0; dex_23=5; dex_31=5;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==1)
		 	{
		 		dex_12=1; dex_23=1; dex_31=1;
		 	}
		 	if (dex_1==2 && dex_2==2 && dex_3==1)
		 	{
		 		dex_12=2; dex_23=3; dex_31=3;
		 	}
		 	if (dex_1==0 && dex_2==0 && dex_3==2)
		 	{
		 		dex_12=0; dex_23=4; dex_31=4;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==2)
		 	{
		 		dex_12=1; dex_23=3; dex_31=3;
		 	}
		 	if (dex_1==2 && dex_2==2 && dex_3==2)
		 	{
		 		dex_12=2; dex_23=2; dex_31=2;
		 	}
		 	break;
		}
		default:
		{
			ExceptionT::BadInputValue(caller, " unsupported spatial dimensions %d", nsd);
		}
	}			
}

//--------------------------------------------------------------------
/* returns jacobian of the nodal values with respect
* to the variables of the shape function second derivatives */
void MFGP_ToolsT::JacobianD2(const LocalArrayT& nodal, const dArray2DT& DDNa,
	dMatrixT& jac) const
{
#if __option(extended_errorcheck)
	/* consistency check */
	const char caller[] = "MFGP_ToolsT::JacobianD2";
	if (nodal.MinorDim() == 1) 
	ExceptionT::OutOfRange(caller, "bad index %d", nodal.MinorDim());
#endif

	double *pjac = jac.Pointer();
	const double *pval = nodal.Pointer();
	
	/* dimensions */
	int nnd   = nodal.NumberOfNodes();
	int nsd = nodal.MinorDim(); //   nodal:[nnd] x [nu] (assumption: nu = nsd)
	int nstr = dSymMatrixT::NumValues(nsd);
	
	/* allocate output space */
	jac.Dimension(nsd, nstr);
	if (nsd == 2)
	{
		double& j11 = *pjac++;
		double& j21 = *pjac++;
		double& j12 = *pjac++;
		double& j22 = *pjac++;
		double& j13 = *pjac++;
		double& j23 = *pjac  ;
		
		j11 = j21 = j12 = j22 = j13 = j23 = 0.0;
		
		const double* pu1 = nodal(0);
		const double* pu2 = nodal(1);
		const double* dxx = DDNa(0);
		const double* dyy = DDNa(1);
		const double* dxy = DDNa(2);
			
		for (int i = 0; i < nnd; i++)
		{
			j11 += (*pu1)*(*dxx);
	   		j21 += (*pu2)*(*dxx);
			j12 += (*pu1)*(*dyy);
			j22 += (*pu2)*(*dyy);
			j13 += (*pu1)*(*dxy);
	   		j23 += (*pu2)*(*dxy);
			
			pu1++; pu2++; dxx++; dyy++; dxy++;	
		}
	}
	else if (nsd == 3)
	{
		double& j11 = *pjac++;
		double& j21 = *pjac++;
		double& j31 = *pjac++;
		double& j12 = *pjac++;
		double& j22 = *pjac++;
		double& j32 = *pjac++;
		double& j13 = *pjac++;
		double& j23 = *pjac++;
		double& j33 = *pjac++;
		double& j14 = *pjac++;
		double& j24 = *pjac++;
		double& j34 = *pjac++;
		double& j15 = *pjac++;
		double& j25 = *pjac++;
		double& j35 = *pjac++;
		double& j16 = *pjac++;
		double& j26 = *pjac++;
		double& j36 = *pjac  ;
	
		j11 = j21 = j31 = j12 = j22 = j32 = j13 = j23 = j33 = 0.0;
		j14 = j24 = j34 = j15 = j25 = j35 = j16 = j26 = j36 = 0.0;

		const double* pu1 = nodal(0);
		const double* pu2 = nodal(1);
		const double* pu3 = nodal(2);
		const double* dxx = DDNa(0);
		const double* dyy = DDNa(1);
		const double* dzz = DDNa(2);
		const double* dyz = DDNa(3);
		const double* dxz = DDNa(4);
		const double* dxy = DDNa(5);

		for (int i = 0; i < nnd; i++)
		{
			j11 += (*pu1)*(*dxx);
			j21 += (*pu2)*(*dxx);
			j31 += (*pu3)*(*dxx);
			j12 += (*pu1)*(*dyy);
			j22 += (*pu2)*(*dyy);
			j32 += (*pu3)*(*dyy);
			j13 += (*pu1)*(*dzz);
			j23 += (*pu2)*(*dzz);
			j33 += (*pu3)*(*dzz);
			j14 += (*pu1)*(*dyz);
			j24 += (*pu2)*(*dyz);
			j34 += (*pu3)*(*dyz);
			j15 += (*pu1)*(*dxz);
			j25 += (*pu2)*(*dxz);
			j35 += (*pu3)*(*dxz);
			j16 += (*pu1)*(*dxy);
			j26 += (*pu2)*(*dxy);
			j36 += (*pu3)*(*dxy);
			
			pu1++; pu2++; pu3++; dxx++; dyy++; dzz++;
			dyz++; dxz++; dxy++;	
		}
	}
	else
	{
		cout << "\n MFGP_ToolsT::JacobianD2: invalid nsd " << endl;
		throw ExceptionT::kBadInputValue;
	}
	
}

//--------------------------------------------------------------------
/* returns jacobian of the nodal values with respect
* to the variables of the shape function third derivatives */
void MFGP_ToolsT::JacobianD3(const LocalArrayT& nodal, const dArray2DT& DDDNa,
	dMatrixT& jac) const
{
#if __option(extended_errorcheck)
	/* consistency check */
	const char caller[] = "MFGP_ToolsT::JacobianD3";
	if (nodal.MinorDim() == 1) 
	ExceptionT::OutOfRange(caller, "bad index %d", nodal.MinorDim());
#endif

	double *pjac = jac.Pointer();
	const double *pval = nodal.Pointer();
	
	/* dimensions */
	int nnd   = nodal.NumberOfNodes();
	int nsd = nodal.MinorDim(); //   nodal:[nnd] x [nu] (assumption: nu = nsd)
	
	jac.Dimension(nsd, nsd*nsd);
	if (nsd == 2)
	{
		double& j11 = *pjac++;
		double& j21 = *pjac++;
		double& j12 = *pjac++;
		double& j22 = *pjac++;
		double& j13 = *pjac++;
		double& j23 = *pjac++;
		double& j14 = *pjac++;
		double& j24 = *pjac;
	
		j11 = j21 = j12 = j22 = 0.0;
		j13 = j23 = j14 = j24 = 0.0;

		const double* pu1 = nodal(0);
		const double* pu2 = nodal(1);
		const double* dxxx = DDDNa(0);
		const double* dyyx = DDDNa(1);
		const double* dxxy = DDDNa(2);
		const double* dyyy = DDDNa(3);

		for (int i = 0; i < nnd; i++)
		{
			j11 += (*pu1)*(*dxxx);
	   		j21 += (*pu2)*(*dxxx);
			j12 += (*pu1)*(*dyyx);
			j22 += (*pu2)*(*dyyx);
			j13 += (*pu1)*(*dxxy);
	   		j23 += (*pu2)*(*dxxy);
			j14 += (*pu1)*(*dyyy);
			j24 += (*pu2)*(*dyyy);
			
			pu1++; pu2++; dxxx++; dyyx++; dxxy++; dyyy++;	
		}
	}
	else if (nsd == 3)
	{
		double& j11 = *pjac++;
		double& j21 = *pjac++;
		double& j31 = *pjac++;
		double& j12 = *pjac++;
		double& j22 = *pjac++;
		double& j32 = *pjac++;
		double& j13 = *pjac++;
		double& j23 = *pjac++;
		double& j33 = *pjac++;
		double& j14 = *pjac++;
		double& j24 = *pjac++;
		double& j34 = *pjac++;
		double& j15 = *pjac++;
		double& j25 = *pjac++;
		double& j35 = *pjac++;
		double& j16 = *pjac++;
		double& j26 = *pjac++;
		double& j36 = *pjac++;
		double& j17 = *pjac++;
		double& j27 = *pjac++;
		double& j37 = *pjac++;
		double& j18 = *pjac++;
		double& j28 = *pjac++;
		double& j38 = *pjac++;
		double& j19 = *pjac++;
		double& j29 = *pjac++;
		double& j39 = *pjac  ;
	
		j11 = j21 = j31 = j12 = j22 = j32 = j13 = j23 = j33 = 0.0;
		j14 = j24 = j34 = j15 = j25 = j35 = j16 = j26 = j36 = 0.0;
		j17 = j27 = j37 = j18 = j28 = j38 = j19 = j29 = j39 = 0.0;

		const double* pu1 = nodal(0);
		const double* pu2 = nodal(1);
		const double* pu3 = nodal(2);
		const double* dxxx = DDDNa(0);
		const double* dyyx = DDDNa(1);
		const double* dzzx = DDDNa(2);
		const double* dxxy = DDDNa(3);
		const double* dyyy = DDDNa(4);
		const double* dzzy = DDDNa(5);
		const double* dxxz = DDDNa(6);
		const double* dyyz = DDDNa(7);
		const double* dzzz = DDDNa(8);

		for (int i = 0; i < nnd; i++)
		{
			j11 += (*pu1)*(*dxxx);
			j21 += (*pu2)*(*dxxx);
			j31 += (*pu3)*(*dxxx);
			j12 += (*pu1)*(*dyyx);
			j22 += (*pu2)*(*dyyx);
			j32 += (*pu3)*(*dyyx);
			j13 += (*pu1)*(*dzzx);
			j23 += (*pu2)*(*dzzx);
			j33 += (*pu3)*(*dzzx);
			j14 += (*pu1)*(*dxxy);
			j24 += (*pu2)*(*dxxy);
			j34 += (*pu3)*(*dxxy);
			j15 += (*pu1)*(*dyyy);
			j25 += (*pu2)*(*dyyy);
			j35 += (*pu3)*(*dyyy);
			j16 += (*pu1)*(*dzzy);
			j26 += (*pu2)*(*dzzy);
			j36 += (*pu3)*(*dzzy);
			j17 += (*pu1)*(*dxxz);
			j27 += (*pu2)*(*dxxz);
			j37 += (*pu3)*(*dxxz);
			j18 += (*pu1)*(*dyyz);
			j28 += (*pu2)*(*dyyz);
			j38 += (*pu3)*(*dyyz);
			j19 += (*pu1)*(*dzzz);
			j29 += (*pu2)*(*dzzz);
			j39 += (*pu3)*(*dzzz);
			
			pu1++; pu2++; pu3++; dxxx++; dyyx++; dzzx++;
			dxxy++; dyyy++; dzzy++; dxxz++; dyyz++; dzzz++;	
		}
	}
	else
	{
		cout << "\n MFGP_ToolsT::JacobianD3: invalid nsd " << endl;
		throw ExceptionT::kBadInputValue;
	}
}