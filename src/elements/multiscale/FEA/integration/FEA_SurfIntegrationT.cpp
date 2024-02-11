// $Id: FEA_SurfIntegrationT.cpp,v 1.2 2003-10-12 23:47:51 raregue Exp $
#include "FEA.h" 

using namespace Tahoe;

FEA_SurfIntegrationT::FEA_SurfIntegrationT() { };

//---------------------------------------------------------

FEA_SurfIntegrationT::FEA_SurfIntegrationT	(FEA_dScalarT &j, FEA_dScalarT &Weights) 
{ 
n_ip 	= Weights.IPs();
W.FEA_Dimension (n_ip);
J.FEA_Dimension (n_ip);
W 		= Weights;
J 		= j; 
};


//---------------------------------------------------------

void FEA_SurfIntegrationT::Construct (FEA_dScalarT &j, FEA_dScalarT &Weights) 
{ 
n_ip 	= Weights.IPs(); 
W.FEA_Dimension (n_ip);
J.FEA_Dimension (n_ip);
W 		= Weights;
J 		= j; 
};

//########################################################## 
//########################################################## 
//########################################################## 


dMatrixT FEA_SurfIntegrationT::of ( FEA_dVectorT &B1, double &c, FEA_dVectorT &B2 )
{
	n_rows = B1.Rows(); // Since using it's transpose
	n_cols = B2.Rows();
	FEA_dMatrixT K(n_ip,n_rows,n_cols); 
	dMatrixT k(n_rows,n_cols); 
	K.Outer(B1,B2);
	K *= c;
	K *= J;
	K *= W;
	
	k = K[0];
	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

	return k;
}



dArrayT FEA_SurfIntegrationT::of ( FEA_dVectorT &B1, double &c, FEA_dScalarT &s)
{
	n_rows = B1.Rows();
	FEA_dVectorT F(n_ip,n_rows); 
	dArrayT f(n_rows); 
	F = B1;
	F *= c;
	F *= s;
	F *= J;
	F *= W;
	
	f = F[0];
	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

	return f;
}



