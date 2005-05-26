/* $Id: MFGP_ToolsT.cpp,v 1.3 2005-05-19 22:05:54 kyonten Exp $ */
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
	const char caller[] = "MFGP_ToolsT::ExpandIndex2";

	switch (nsd)
	{
		case 1:
		{
		 if (dex_1==0 && dex_2==0 && dex_3==0){
		 	dex_12=0; dex_23=0; dex_31=0;
		 }
		 break;	
		}
		
		case 2:
		{
			if (dex_1==0 && dex_2==0 && dex_3==0){
		 		dex_12=0; dex_23=0; dex_31=0;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==0){
		 		dex_12=1; dex_23=2; dex_31=2;
		 	}
		 	if (dex_1==0 && dex_2==0 && dex_3==1){
		 		dex_12=0; dex_23=2; dex_31=2;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==1){
		 		dex_12=1; dex_23=1; dex_31=1;
		 	}
		 	break;
		}
		
		case 3:
		{
			if (dex_1==0 && dex_2==0 && dex_3==0){
		 		dex_12=0; dex_23=0; dex_31=0;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==0){
		 		dex_12=1; dex_23=5; dex_31=5;
		 	}
		 	if (dex_1==2 && dex_2==2 && dex_3==0){
		 		dex_12=2; dex_23=4; dex_31=4;
		 	}
		 	if (dex_1==0 && dex_2==0 && dex_3==1){
		 		dex_12=0; dex_23=5; dex_31=5;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==1){
		 		dex_12=1; dex_23=1; dex_31=1;
		 	}
		 	if (dex_1==2 && dex_2==2 && dex_3==1){
		 		dex_12=2; dex_23=3; dex_31=3;
		 	}
		 	if (dex_1==0 && dex_2==0 && dex_3==2){
		 		dex_12=0; dex_23=4; dex_31=4;
		 	}
		 	if (dex_1==1 && dex_2==1 && dex_3==2){
		 		dex_12=1; dex_23=3; dex_31=3;
		 	}
		 	if (dex_1==2 && dex_2==2 && dex_3==2){
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
