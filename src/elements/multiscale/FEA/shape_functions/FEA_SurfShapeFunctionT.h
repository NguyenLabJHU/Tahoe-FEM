// $Id: FEA_SurfShapeFunctionT.h,v 1.1 2003-10-07 19:37:49 raregue Exp $
#ifndef _FEA_SURFSHAPEFUNCTIONT_H_
#define _FEA_SURFSHAPEFUNCTIONT_H_

#include "FEA.h"

namespace Tahoe {

//------------------------------------------------------------------------
				
class FEA_SurfShapeFunctionT
{
	public:
	
    	FEA_SurfShapeFunctionT 	() { n_ip=n_sd=n_en=0; }
    	FEA_SurfShapeFunctionT 	(int &nip,int &nsd,int &nen); 
    	void Construct 		(int &nip,int &nsd,int &nen); 
		void Print 			( char* name);

		FEA_dVectorT N;						
		FEA_dMatrixT dNdx;						
		FEA_dScalarT W;						
		FEA_dScalarT j;						

		int n_ip, n_sd, n_en;
};

//------------------------------------------------------------------------

inline FEA_SurfShapeFunctionT::FEA_SurfShapeFunctionT (int &nip,int &nsd,int &nen) 
{ 
	n_ip=nip; n_sd=nsd; n_en=nen; 

	     N.FEA_Dimension	(	nip,nen	); 
  		dNdx.FEA_Dimension	(	nip,nsd,nen	); 
 	     W.FEA_Dimension	(	nip	); 
 	     j.FEA_Dimension	(	nip	); 
}

//------------------------------------------------------------------------

inline void FEA_SurfShapeFunctionT::Construct (int &nip,int &nsd,int &nen) 
{ 
	n_ip=nip; n_sd=nsd; n_en=nen; 

	     N.FEA_Dimension	(	nip,nen	); 
  		dNdx.FEA_Dimension	(	nip,nsd,nen	); 
 	     W.FEA_Dimension	(	nip	); 
 	     j.FEA_Dimension	(	nip	); 
}

//------------------------------------------------------------------------

inline void FEA_SurfShapeFunctionT::Print (char* name)
{ 
	cout << " FEA_SurfShapeFunctionT " << name << "\n\n";
	N.Print("N");
	dNdx.Print("dNdx");
	W.Print("W");
	j.Print("j");

}

}

#endif
