// $Id: APS_FEA_Data_Processor_DisplT.cpp,v 1.1 2004-02-17 19:48:49 raregue Exp $
#include "APS_FEA_Data_Processor_DisplT.h"  

using namespace Tahoe;

//---------------------------------------------------------------------

APS_FEA_Data_Processor_DisplT::APS_FEA_Data_Processor_DisplT() { };

//---------------------------------------------------------------------

APS_FEA_Data_Processor_DisplT::~APS_FEA_Data_Processor_DisplT() { };

//---------------------------------------------------------------------

APS_FEA_Data_Processor_DisplT::APS_FEA_Data_Processor_DisplT( FEA_dMatrixT &fdNdx ) 
{
	Construct ( fdNdx );
}

//---------------------------------------------------------------------

void APS_FEA_Data_Processor_DisplT::Construct ( FEA_dMatrixT &fdNdx )  
{
	if (fdNdx.Length()==0)  
 		cout <<"...ERROR >> APS_FEA_Data_Processor_DisplT::Initialize_Data_Pro : dNadx unallocated \n\n";

	dN = fdNdx;
  	n_ip = dN.IPs();
  	n_sd = dN.Rows(); 
  	n_en = dN.Cols();
}

//---------------------------------------------------------------------

void APS_FEA_Data_Processor_DisplT::APS_B	(FEA_dMatrixT &B) 
{
			B = 0.0;
		
			for (int a=0; a<n_en; a++)
	  			B(0,a) = dN(dx1,a); 
	  						
	  		for (int a=0; a<n_en; a++)
	  			B(1,a) = dN(dx2,a); 	  			
}


