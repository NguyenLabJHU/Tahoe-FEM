// $Id: APS_FEA_Data_Processor_SurfT.cpp,v 1.2 2003-10-08 17:45:16 raregue Exp $
#include "APS_FEA_Data_Processor_SurfT.h"  

using namespace Tahoe;

//---------------------------------------------------------------------

APS_FEA_Data_Processor_SurfT::APS_FEA_Data_Processor_SurfT() { };

//---------------------------------------------------------------------

APS_FEA_Data_Processor_SurfT::~APS_FEA_Data_Processor_SurfT() { };

//---------------------------------------------------------------------

APS_FEA_Data_Processor_SurfT::APS_FEA_Data_Processor_SurfT( FEA_dMatrixT &fdNdx_surf ) 
{
	Construct ( fdNdx_surf );
}

//---------------------------------------------------------------------

void APS_FEA_Data_Processor_SurfT::Construct ( FEA_dMatrixT &fdNdx_surf )  
{
  if (fdNdx_surf.Length()==0)  
 		cout <<"...ERROR >> APS_FEA_Data_Processor_SurfT::Initialize_Data_Pro : dNadx unallocated \n\n";

	dN_surf = fdNdx_surf;
  	n_ip = dN_surf.IPs();
  	n_sd = dN_surf.Rows(); 
  	n_en = dN_surf.Cols();
}

//---------------------------------------------------------------------

void APS_FEA_Data_Processor_SurfT::APS_B_surf (FEA_dMatrixT &B) 
{
	B = 0.0;
	int a;
	for (a=0; a<n_en; a++) B(0,a) = dN_surf(dx1,a); 
	for (a=0; a<n_en; a++) B(1,a) = dN_surf(dx2,a); 
}

//---------------------------------------------------------------------

void APS_FEA_Data_Processor_SurfT::APS_N	(FEA_dVectorT &B) 
{
	B = 0.0;
	int a;
	for (a=0; a<n_en; a++) B(a) = N_surf(a); 	
}

