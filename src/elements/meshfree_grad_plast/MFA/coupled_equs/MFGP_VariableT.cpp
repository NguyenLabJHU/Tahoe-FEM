// $Id: MFGP_VariableT.cpp
#include "MFGP_VariableT.h"

//---------------------------------------------------------------------
/** constructor */

using namespace Tahoe;

MFGP_VariableT::MFGP_VariableT (const FEA_dMatrixT& grad_u, const FEA_dVectorT& gammap, 
							  const FEA_dMatrixT& grad_gammap, FEA_dVectorT& state)
{
	Construct	(grad_u, gammap, grad_gammap, state);
}

//---------------------------------------------------------------------
/** Data initialization: Allocate space for grad_u, gammap, grad_gammap  */
 
void MFGP_VariableT::Construct (const FEA_dMatrixT& grad_u, const FEA_dVectorT& gammap,
							   const FEA_dMatrixT& grad_gammap, FEA_dVectorT& state)
{
  	n_vars_vector = MFGP::kNUM_MFGP_VECTOR_VARS;  
  	n_vars_matrix = MFGP::kNUM_MFGP_MATRIX_VARS; 
  	fVars_vector.Dimension( n_vars_vector );
  	fVars_matrix.Dimension( n_vars_matrix );  

	fVars_matrix[MFGP::kgrad_u] = grad_u;  
	fVars_vector[MFGP::kgammap] = gammap;  
	fVars_matrix[MFGP::kgrad_gammap] = grad_gammap; 
	fVars_vector[MFGP::kstate] = state;
}

//----------------------------------------------------

void MFGP_VariableT::Delete_Vars	( void )
{
#if 0
	for (int i=0; i<n_vars; i++) 
		fVars_vector[i].FEA_Delete(); // ArrayT checks if fLength=0 before deleting
#endif
#pragma message("MFGP_VariableT::Delete_Vars: fVars[i] is a FEA_dVectorT which has no FEA_Delete()")
}

//----------------------------------------------------

void MFGP_VariableT::Print() { Print(" "); } 

void MFGP_VariableT::Print(char *c) { // overload << later

  cout <<"\n MFGP_VariableT "<< c <<" follows: \n\n"; 

	for (int l=0; l<n_vars_vector; l++) 
		if (fVars_vector[l].IPs() == 0)
			cout << "MFGP_VariableT n["<<l<<"] Unallocated \n\n";
    else {
  		cout << "\n Vector "<<l<<" evaluated at "<<fVars_vector[l].IPs() <<" integration points (ip): \n"; 
  		for (int i=0; i<fVars_vector[l].IPs(); i++) 
   	 		cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< fVars_vector[l][i] << "\n";
			cout << "\n";
		}
		
	for (int l=0; l<n_vars_matrix; l++) 
		if (fVars_matrix[l].IPs() == 0)
			cout << "MFGP_VariableT n["<<l<<"] Unallocated \n\n";
    else {
  		cout << "\n Vector "<<l<<" evaluated at "<<fVars_matrix[l].IPs() <<" integration points (ip): \n"; 
  		for (int i=0; i<fVars_matrix[l].IPs(); i++) 
   	 		cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< fVars_matrix[l][i] << "\n";
			cout << "\n";
		}
	
}

//---------------------------------------------------------------------
//** Put variable from local model calculation to class work space
void MFGP_VariableT::Put(MFGP::VarT_vector variable, FEA_dVectorT& var)  
{
  fVars_vector[variable] = var;
} 

//---------------------------------------------------------------------
//** Update variable from class work space
void MFGP_VariableT::Update(MFGP::VarT_vector variable, FEA_dVectorT& var)  
{
  var = fVars_vector[variable];
} 

//---------------------------------------------------------------------
//** Retrieve/Fetch/Get either grad_u or gammap from class work space
const FEA_dVectorT& MFGP_VariableT::Get(MFGP::VarT_vector variable)  
{
  if (!fVars_vector[variable].Length())  // 0 rows indicated un-allocation
    Allocate_and_Compute_Variables(variable);

  return fVars_vector[variable];
} 

//---------------------------------------------------------------------
//** Retrieve/Fetch/Get either grad_u or gammap from class work space
const FEA_dMatrixT& MFGP_VariableT::Get(MFGP::VarT_matrix variable)  
{
  if (!fVars_matrix[variable].Length())  // 0 rows indicated un-allocation
    Allocate_and_Compute_Variables(variable);

  return fVars_matrix[variable];
} 


//---------------------------------------------------------------------
/** Compute and store ... : recursive routine */ 
void MFGP_VariableT::Allocate_and_Compute_Variables(MFGP::VarT_vector kVariable)
{
    switch (kVariable) 
    {

	case MFGP::kgammap : // gammap   
        fVars_vector[MFGP::kgammap]; 
				break;
				
	case MFGP::kstate : // ISVs   
        fVars_vector[MFGP::kstate]; 
				break;		
				
	default:
        cout << "\n MFGP_VariableT:::Allocate_and_Compute_Variables() bad VarT type" << endl;
				//throw eGeneralFail;
				break;
    }
}


//---------------------------------------------------------------------
/** Compute and store ... : recursive routine */ 
void MFGP_VariableT::Allocate_and_Compute_Variables(MFGP::VarT_matrix kVariable)
{
    switch (kVariable) 
    {
    
	case MFGP::kgrad_u : // grad_u 
        fVars_matrix[MFGP::kgrad_u]; 
				break;
	
				
	case MFGP::kgrad_gammap : // grad_gammap   
        fVars_matrix[MFGP::kgrad_gammap]; 
				break;

	default:
        cout << "\n MFGP_VariableT:::Allocate_and_Compute_Variables() bad VarT type" << endl;
				//throw eGeneralFail;
				break;
    }
}


//---------------------------------------------------------------------

void MFGP_VariableT::operator=(const MFGP_VariableT &a)	// Initializes
{
	n_vars_vector = a.n_vars_vector;
	fVars_vector[ n_vars_vector ];

	for (int i=0; i<n_vars_vector; i++) 
	{
  	fVars_vector[i].FEA_Dimension( a.fVars_vector[i].IPs(), a.fVars_vector[i].Rows() );
  	fVars_vector[i] = a.fVars_vector[i];
	}

	n_vars_matrix = a.n_vars_matrix;
	fVars_matrix.Dimension( n_vars_matrix );

	for (int i=0; i<n_vars_matrix; i++) 
	{
  	fVars_matrix[i].FEA_Dimension( a.fVars_matrix[i].IPs(), a.fVars_matrix[i].Rows(), a.fVars_matrix[i].Cols() );
  	fVars_matrix[i] = a.fVars_matrix[i];
	}
};

//---------------------------------------------------------------------

void MFGP_VariableT::operator +=  (const double &a)
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] += a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] += a;
};


//---------------------------------------------------------------------

void MFGP_VariableT::operator -=  (const double &a)
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] -= a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] -= a;
};

//---------------------------------------------------------------------

void MFGP_VariableT::operator *=  (const double &a)
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] *= a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] *= a;
};

//---------------------------------------------------------------------

void MFGP_VariableT::operator /=  (const double &a) 
{
for (int i=0; i<n_vars_vector; i++) fVars_vector[i] /= a;
for (int i=0; i<n_vars_matrix; i++) fVars_matrix[i] /= a;
};

//---------------------------------------------------------------------

