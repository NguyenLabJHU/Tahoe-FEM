// $Id: APS_VariableT.cpp,v 1.9 2003-09-25 20:40:20 raregue Exp $
#include "APS_VariableT.h"

//---------------------------------------------------------------------
/** constructor */

using namespace Tahoe;

//APS_VariableT::APS_VariableT (const FEA_dVectorT& grad_u, const FEA_dVectorT& gammap, const FEA_dMatrixT& grad_gammap) 
APS_VariableT::APS_VariableT (const FEA_dMatrixT& grad_u, const FEA_dVectorT& gammap, const FEA_dMatrixT& grad_gammap)
{
	Construct	(grad_u, gammap, grad_gammap);
}

//---------------------------------------------------------------------
/** Data initialization: Allocate space for grad_u, gammap, grad_gammap  */

//void APS_VariableT::Construct (const FEA_dVectorT& grad_u, const FEA_dVectorT& gammap,const FEA_dMatrixT& grad_gammap) 
void APS_VariableT::Construct (const FEA_dMatrixT& grad_u, const FEA_dVectorT& gammap,const FEA_dMatrixT& grad_gammap)
{
  	n_vars_vector = APS::kNUM_APS_VECTOR_VARS;  
  	n_vars_matrix = APS::kNUM_APS_MATRIX_VARS; 
  	fVars_vector.Dimension( n_vars_vector );
  	fVars_matrix.Dimension( n_vars_matrix );  

	//fVars_vector[APS::kgrad_u] = grad_u; // This = opr allocates if LHS Length=0
	fVars_matrix[APS::kgrad_u] = grad_u; // This = opr allocates if LHS Length=0
	fVars_vector[APS::kgammap] = gammap; // This = opr allocates if LHS Length=0
	fVars_matrix[APS::kgrad_gammap] = grad_gammap; 
}

//----------------------------------------------------

void APS_VariableT::Delete_Vars	( void )
{
#if 0
	for (int i=0; i<n_vars; i++) 
		fVars_vector[i].FEA_Delete(); // ArrayT checks if fLength=0 before deleting
#endif
#pragma message("APS_VariableT::Delete_Vars: fVars[i] is a FEA_dVectorT which has no FEA_Delete()")
}

//----------------------------------------------------

void APS_VariableT::Print() { Print(" "); } 

void APS_VariableT::Print(char *c) { // overload << later

  cout <<"\n APS_VariableT "<< c <<" follows: \n\n"; 

	for (int l=0; l<n_vars_vector; l++) 
		if (fVars_vector[l].IPs() == 0)
			cout << "APS_VariableT n["<<l<<"] Unallocated \n\n";
    else {
  		cout << "\n Vector "<<l<<" evaluated at "<<fVars_vector[l].IPs() <<" integration points (ip): \n"; 
  		for (int i=0; i<fVars_vector[l].IPs(); i++) 
   	 		cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< fVars_vector[l][i] << "\n";
			cout << "\n";
		}
		
	for (int l=0; l<n_vars_matrix; l++) 
		if (fVars_matrix[l].IPs() == 0)
			cout << "APS_VariableT n["<<l<<"] Unallocated \n\n";
    else {
  		cout << "\n Vector "<<l<<" evaluated at "<<fVars_matrix[l].IPs() <<" integration points (ip): \n"; 
  		for (int i=0; i<fVars_matrix[l].IPs(); i++) 
   	 		cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< fVars_matrix[l][i] << "\n";
			cout << "\n";
		}
	
}

//---------------------------------------------------------------------
//** Retrieve/Fetch/Get either grad_u or gammap from class work space
const FEA_dVectorT& APS_VariableT::Get(APS::VarT_vector variable)  
{

  if (!fVars_vector[variable].Length())  // 0 rows indicated un-allocation
    Allocate_and_Compute_Variables(variable);

  return fVars_vector[variable];
  
} 

//---------------------------------------------------------------------
//** Retrieve/Fetch/Get either grad_u or gammap from class work space
const FEA_dMatrixT& APS_VariableT::Get(APS::VarT_matrix variable)  
{

  if (!fVars_matrix[variable].Length())  // 0 rows indicated un-allocation
    Allocate_and_Compute_Variables(variable);

  return fVars_matrix[variable];
  
} 


//---------------------------------------------------------------------
/** Compute and store ... : recursive routine */ 
void APS_VariableT::Allocate_and_Compute_Variables(APS::VarT_vector kVariable)
{

    switch (kVariable) {

      /*case APS::kgrad_u : // grad_u 
        fVars_vector[APS::kgrad_u]; 
				break;
				*/

      case APS::kgammap : // gammap   
        fVars_vector[APS::kgammap]; 
				break;
				
      default:
        cout << "\n APS_VariableT:::Allocate_and_Compute_Variables() bad VarT type" << endl;
				//throw eGeneralFail;
				break;
    }
}


//---------------------------------------------------------------------
/** Compute and store ... : recursive routine */ 
void APS_VariableT::Allocate_and_Compute_Variables(APS::VarT_matrix kVariable)
{

    switch (kVariable) {
    
    	case APS::kgrad_u : // grad_u 
        fVars_matrix[APS::kgrad_u]; 
				break;
				
      case APS::kgrad_gammap : // grad_gammap   
        fVars_matrix[APS::kgrad_gammap]; 
				break;

      default:
        cout << "\n APS_VariableT:::Allocate_and_Compute_Variables() bad VarT type" << endl;
				//throw eGeneralFail;
				break;
    }
}


//---------------------------------------------------------------------

void APS_VariableT::operator=(const APS_VariableT &a)	// Initializes
{

	n_vars_vector = a.n_vars_vector;
	fVars_vector[ n_vars_vector ];

	for (int i=0; i<n_vars_vector; i++) {
  	fVars_vector[i].FEA_Dimension( a.fVars_vector[i].IPs(), a.fVars_vector[i].Rows() );
  	fVars_vector[i] = a.fVars_vector[i];
	}

	n_vars_matrix = a.n_vars_matrix;
	fVars_matrix.Dimension( n_vars_matrix );

	for (int i=0; i<n_vars_matrix; i++) {
  	fVars_matrix[i].FEA_Dimension( a.fVars_matrix[i].IPs(), a.fVars_matrix[i].Rows(), a.fVars_matrix[i].Cols() );
  	fVars_matrix[i] = a.fVars_matrix[i];
	}
};

//---------------------------------------------------------------------

void APS_VariableT::operator +=  (const double &a)
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] += a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] += a;
};


//---------------------------------------------------------------------

void APS_VariableT::operator -=  (const double &a)
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] -= a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] -= a;
};

//---------------------------------------------------------------------

void APS_VariableT::operator *=  (const double &a)
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] *= a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] *= a;
};

//---------------------------------------------------------------------

void APS_VariableT::operator /=  (const double &a) 
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] /= a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] /= a;
};

//---------------------------------------------------------------------

