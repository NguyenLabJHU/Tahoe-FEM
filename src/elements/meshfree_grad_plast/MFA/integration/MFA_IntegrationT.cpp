// $Id: MFA_IntegrationT.cpp
#include "MFA.h" //This class doesn't exist. see multiscale folder

using namespace Tahoe;

MFA_IntegrationT::MFA_IntegrationT() { };

//---------------------------------------------------------

MFA_IntegrationT::MFA_IntegrationT	(double Weights) 
{ 
n_ip 	= Weights.IPs();
W.Dimension (n_ip);
W 		= Weights;
};

//---------------------------------------------------------


//########################################################## 
//########################################################## 

// LHS integrals
//########################################################## 

dMatrixT MFA_IntegrationT::of ( dMatrixT &B1, dMatrixT &C, dMatrixT &B2 )
{
  n_rows = B1.Cols(); // Since using its transpose
  n_cols = B2.Cols();
	dMatrixT 	K	(n_ip,n_rows,n_cols); 
	dMatrixT 			k(n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------
dArrayT MFA_IntegrationT::of ( dMatrixT &B1, dArrayT &C_v, double c )
{
  n_rows = B1.Cols(); // since using its transpose
	dArrayT 	F	(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B1,C_v);
	F *= c;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------
dArrayT MFA_IntegrationT::of ( double c, dArrayT &C_v, dMatrixT &B1)
{
  n_rows = B1.Cols(); // since using its transpose
	dArrayT 	F	(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B1,C_v);
	F *= c;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}


//---------------------------------------------------------
double MFA_IntegrationT::of ( double c, double cc, double ccc )
{
	int constant = 1;
	dArrayT 	K(n_ip,constant); //not too sure, double check!!
	double k;
	K = c * cc;
	K *= ccc;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}


//RHS integrals
//########################################################## 

dArrayT MFA_IntegrationT::of ( dMatrixT &B, dArrayT &b )
{
  n_rows = B.Cols(); // Since using its transpose
	dArrayT 	F(n_ip,n_rows); 
	dArrayT 			f	(n_rows); 
	F.MultATb 			(B,b);
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------

dArrayT MFA_IntegrationT::of ( double c, dArrayT &b )
{
  n_rows = b.Cols(); // Since using its transpose
	dArrayT 	F(n_ip,n_rows); 
	dArrayT 			f(n_rows); 
	F = c * b;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}


//---------------------------------------------------------

double MFA_IntegrationT::of ( double c, double cc)
{
	int constant = 1;
	dArrayT 	K	(n_ip,constant); //not sure, double check!!
	double k;
	K = c * cc;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}


//########################################################################################################### 
//########################################################################################################### 

//LHS integrals
//########################################################## 

void MFA_IntegrationT::of ( dMatrixT &B1, dMatrixT &C, dMatrixT &B2)
{
  n_rows = B1.Cols(); // Since using it's transpose
  n_cols = B2.Cols();
	dMatrixT 	K	(n_ip,n_rows,n_cols); 
	K.MultATBC 			(B1,C,B2);
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

}
  	
//---------------------------------------------------------
void MFA_IntegrationT::of ( dMatrixT &B1, dArrayT &C_v, double c )
{
  n_rows = B1.Cols(); // since using its transpose
	dArrayT 	F(n_ip,n_rows);  
	F.MultATb(B1,C_v);
	F *= c;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}


//---------------------------------------------------------
void MFA_IntegrationT::of ( double c, dArrayT &C_v, dMatrixT &B1)
{
  n_rows = B1.Cols(); // since using its transpose
	FEA_dVectorT 	F	(n_ip,n_rows);  
	F.MultATb 			(B1,C_v);
	F *= c;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}


//---------------------------------------------------------
void MFA_IntegrationT::of ( double c, double cc, double ccc )
{
	int constant = 1;
	dArrayT 	K	(n_ip,constant); //not sure, double check!!
	double k;
	K = c * cc;
	K *= ccc;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}


//RHS integrals
//########################################################## 

void MFA_IntegrationT::of ( dMatrixT &B, dArrayT &b )
{
  n_rows = B.Cols(); // Since using its transpose
	dArrayT 	F	(n_ip,n_rows);  
	F.MultATb 			(B,b);
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}

//---------------------------------------------------------

void MFA_IntegrationT::of ( double c, dArrayT &b )
{
  n_rows = b.Cols(); // Since using its transpose
	dArrayT 	F	(n_ip,n_rows); 
	F = c * b;
	F *= W;
	f = F[0];

	for (int i=1; i<n_ip; i++) // Summation here
		f += F[i];

return f;
}


//---------------------------------------------------------

void MFA_IntegrationT::of ( double c, double cc)
{
	int constant = 1;
	dArrayT 	K	(n_ip,constant); //not sure, double check!!
	double k;
	K = c * cc;
	K *= W;
	k = K[0];

	for (int i=1; i<n_ip; i++) // Summation here
		k += K[i];

return k;
}

//---------------------------------------------------------