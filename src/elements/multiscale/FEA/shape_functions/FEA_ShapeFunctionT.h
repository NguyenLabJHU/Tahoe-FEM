//DEVELOPMENT

#ifndef _FEA_SHAPEFUNCTIONT_H_
#define _FEA_SHAPEFUNCTIONT_H_

#include "FEA.h"

namespace Tahoe {


class FEA_ShapeFunctionT
{
	public:
	
    FEA_ShapeFunctionT () { n_ip=n_sd=n_en=0; }
    FEA_ShapeFunctionT (int &nip,int &nsd,int &nen); 

		FEA_dVectorT N;						
		FEA_dMatrixT dNdx;						
		FEA_dScalarT W;						
		FEA_dScalarT j;						

		int n_ip, n_sd, n_en;
};

//------------------------------------------------------------------------

inline  FEA_ShapeFunctionT::FEA_ShapeFunctionT (int &nip,int &nsd,int &nen) 
{ 
	n_ip=nip; n_sd=nsd; n_en=nen; 

	     N.FEA_Dimension(nip,nen); 
  	dNdx.FEA_Dimension(nip,nsd,nen); 
 	     W.FEA_Dimension(nip); 
 	     j.FEA_Dimension(nip); 
}

}
// inline routines go here ...

#endif
