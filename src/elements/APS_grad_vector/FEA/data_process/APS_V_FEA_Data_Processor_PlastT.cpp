// $Id: APS_V_FEA_Data_Processor_PlastT.cpp,v 1.1 2004-03-02 23:46:52 raregue Exp $
#include "APS_V_FEA_Data_Processor_PlastT.h"  

using namespace Tahoe;

//---------------------------------------------------------------------

APS_V_FEA_Data_Processor_PlastT::APS_V_FEA_Data_Processor_PlastT() { };

//---------------------------------------------------------------------

APS_V_FEA_Data_Processor_PlastT::~APS_V_FEA_Data_Processor_PlastT() { };

//---------------------------------------------------------------------

APS_V_FEA_Data_Processor_PlastT::APS_V_FEA_Data_Processor_PlastT( FEA_dMatrixT &fdNdx ) 
{
	Construct ( fdNdx );
}

//---------------------------------------------------------------------

void APS_V_FEA_Data_Processor_PlastT::Construct ( FEA_dMatrixT &fdNdx )  
{
	if (fdNdx.Length()==0)  
 		cout <<"...ERROR >> APS_V_FEA_Data_Processor_PlastT::Initialize_Data_Pro : dNadx unallocated \n\n";

	dN = fdNdx;
  	n_ip = dN.IPs();
  	n_sd = dN.Rows(); 
  	n_en = dN.Cols();
}


//---------------------------------------------------------------------

void APS_V_FEA_Data_Processor_PlastT::APS_Ngamma (FEA_dMatrixT &B) 
{
			B = 0.0;
			
			for (int a=0; a<n_en; a++)
	  			B(0,a*2) = N(a); 
	  						
	  		for (int a=0; a<n_en; a++)
	  			B(1,(a*2)+1) = N(a); 	  			
}

//---------------------------------------------------------------------

void APS_V_FEA_Data_Processor_PlastT::APS_Ngam1d2	(FEA_dVectorT &B) 
{
			B = 0.0;
			
			for (int a=0; a<n_en; a++) 
	  			B(a*2) = dN(dx2,a);   			
}

//---------------------------------------------------------------------

void APS_V_FEA_Data_Processor_PlastT::APS_Ngam2d1	(FEA_dVectorT &B) 
{
			B = 0.0;
			
	  		for (int a=0; a<n_en; a++)
	  			B((a*2)+1) = dN(dx1,a);	
}
