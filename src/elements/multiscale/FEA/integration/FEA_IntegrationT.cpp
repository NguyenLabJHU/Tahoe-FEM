// $Id: FEA_IntegrationT.cpp,v 1.5 2003-10-06 18:31:52 raregue Exp $
#include "FEA.h" 

using namespace Tahoe;

FEA_IntegrationT::FEA_IntegrationT() { };

//---------------------------------------------------------

FEA_IntegrationT::FEA_IntegrationT	(FEA_dScalarT &j, FEA_dScalarT &Weights) 
{ 
n_ip 	= Weights.IPs();
W.FEA_Dimension (n_ip);
J.FEA_Dimension (n_ip);
W 		= Weights;
J 		= j; 
};


//---------------------------------------------------------

void FEA_IntegrationT::Construct (FEA_dScalarT &j, FEA_dScalarT &Weights) 
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

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &K ) 
{
	dMatrixT k(K.Rows(),K.Cols()); 
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}
  	
//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( double c, FEA_dMatrixT &K )
{
	dMatrixT k(K.Rows(),K.Cols()); 
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dScalarT &s, FEA_dMatrixT &K )
{
	dMatrixT k(K.Rows(),K.Cols()); 
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( double c,FEA_dScalarT &s, FEA_dMatrixT &K )
{
	dMatrixT k(K.Rows(),K.Cols()); 
	K *= c;
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//########################################################## 


dMatrixT FEA_IntegrationT::of ( FEA_dVectorT &B1, double c, FEA_dVectorT &B2 )
{
  n_rows = B1.Rows(); // Since using it's transpose
  n_cols = B2.Rows();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.Outer 			(B1,B2);
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------


dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}
  	
//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, double c, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dScalarT &s, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1,double c,FEA_dScalarT &s, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= c;
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//########################################################## 

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dMatrixT &C, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}
  	
//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, double c, FEA_dMatrixT &C, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dScalarT &s, FEA_dMatrixT &C, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------

dMatrixT FEA_IntegrationT::of ( FEA_dMatrixT &B1, double c, FEA_dScalarT &s, FEA_dMatrixT &C, FEA_dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k	(n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= c;
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//########################################################## 

dArrayT FEA_IntegrationT::of ( FEA_dVectorT &B1, double c)
{
	n_rows = B1.Rows();
	FEA_dVectorT 	F	(n_ip,n_rows); 
	dArrayT 		f	(n_rows); 
	F = B1;
	F *= c;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------

dArrayT FEA_IntegrationT::of ( FEA_dMatrixT &B, FEA_dVectorT &b )
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B,b);
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------

dArrayT FEA_IntegrationT::of ( FEA_dMatrixT &B, double &c, FEA_dVectorT &b )
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B,b);
	F *= c;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------

dArrayT FEA_IntegrationT::of ( FEA_dMatrixT &B, FEA_dScalarT &s, FEA_dVectorT &b )
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B,b);
	F *= s;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------

dArrayT FEA_IntegrationT::of ( FEA_dMatrixT &B,double &c,FEA_dScalarT &s,FEA_dVectorT &b )
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B,b);
	F *= c;
	F *= s;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//########################################################################################################### 
//########################################################################################################### 
//########################################################################################################### 

void FEA_IntegrationT::of ( FEA_dMatrixT &K, dMatrixT &k )
{
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}
  	
//---------------------------------------------------------

void FEA_IntegrationT::of ( double c, FEA_dMatrixT &K, dMatrixT &k )
{
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dScalarT &s, FEA_dMatrixT &K, dMatrixT &k )
{
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//---------------------------------------------------------

void FEA_IntegrationT::of ( double c,FEA_dScalarT &s, FEA_dMatrixT &K, dMatrixT &k )
{
	K *= c;
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//########################################################## 

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dMatrixT &B2, dMatrixT &k )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}
  	
//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, double c, FEA_dMatrixT &B2, dMatrixT &k )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dScalarT &s, FEA_dMatrixT &B2, dMatrixT &k )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B1,double c,FEA_dScalarT &s, FEA_dMatrixT &B2, dMatrixT &k )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATB 			(B1,B2);
	K *= c;
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//########################################################## 

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k )
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}
  	
//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, double c, FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k ) //Untested
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= c;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, FEA_dScalarT &s, FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k ) //Untested
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B1, double c, FEA_dScalarT &s, FEA_dMatrixT &C, FEA_dMatrixT &B2, dMatrixT &k ) //Untested
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	FEA_dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= c;
	K *= s;
	K *= J;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}

//########################################################## 

void FEA_IntegrationT::of ( FEA_dMatrixT &B, FEA_dVectorT &b, dArrayT &f )
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	F.MultATb 			(B,b);
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];
}


//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B, double &c, FEA_dVectorT &b, dArrayT &f ) //Untested
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	F.MultATb 			(B,b);
	F *= c;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];
}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B, FEA_dScalarT &s, FEA_dVectorT &b, dArrayT &f ) //Untested
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	F.MultATb 			(B,b);
	F *= s;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];
}

//---------------------------------------------------------

void FEA_IntegrationT::of ( FEA_dMatrixT &B,double &c,FEA_dScalarT &s,FEA_dVectorT &b, dArrayT &f ) //Untested
{
  n_rows = B.Cols(); // Since using it's transpose
	FEA_dVectorT 	F	(n_ip,n_rows); 
	F.MultATb 			(B,b);
	F *= c;
	F *= s;
	F *= J;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];
}

