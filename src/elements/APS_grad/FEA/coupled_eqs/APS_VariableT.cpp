// $Id: APS_VariableT.cpp,v 1.2 2003-09-02 23:57:54 raregue Exp $
#include "APS_VariableT.h"

//---------------------------------------------------------------------
/** constructor */

using namespace Tahoe;

APS_VariableT::APS_VariableT (const FEA_dVectorT& grad_u, const FEA_dVectorT& gammap) 
{
	Construct	(grad_u, gammap);
}

//---------------------------------------------------------------------
/** Data initialization: Allocate space for grad_u, gammap  */

void APS_VariableT::Construct (const FEA_dVectorT& grad_u, const FEA_dVectorT& gammap) 
{
  	n_vars = APS::kNUM_APS_VARS;  
  	fVars.Dimension( n_vars );  

	fVars[APS::kgrad_u] = grad_u; // This = opr allocates if LHS Length=0
	fVars[APS::kgammap] = gammap; 
}

//----------------------------------------------------

void APS_VariableT::Delete_Vars	( void )
{
	for (int i=0; i<n_vars; i++) 
		fVars[i].FEA_Delete(); // ArrayT checks if fLength=0 before deleting
}

//----------------------------------------------------

void APS_VariableT::Print() { Print(" "); } 

void APS_VariableT::Print(char *c) { // overload << later

  cout <<"\n APS_VariableT "<< c <<" follows: \n\n"; 

	for (int l=0; l<n_vars; l++) 
		if (fVars[l].n_ip==0)
			cout << "APS_VariableT n["<<l<<"] Unallocated \n\n";
    else {
  		cout << "\n Vector "<<l<<" evaluated at "<<fVars[l].n_ip<<" inegration points (ip): \n"; 
  		for (int i=0; i<fVars[l].n_ip; i++) 
   	 		cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< fVars[l][i] << "\n";
			cout << "\n";
		}
	
}

//---------------------------------------------------------------------
//** Retrieve/Fetch/Get either grad_u or gammap from class work space
/*const FEA_dMatrixT& APS_VariableT::Get(APS::VarT variable)  
{

  if (!fVars[variable].Length())  // 0 rows indicated un-allocation
    Allocate_and_Compute_Variables(variable);

  return fVars[variable];
  
}*/

//---------------------------------------------------------------------

void APS_VariableT::operator  =  (const APS_VariableT &a)	// Initializes
{
	n_vars = a.n_vars;
	fVars.Dimension( n_vars );

	for (int i=0; i<n_vars; i++) {
  	fVars[i].FEA_Dimension( a.fVars[i].n_ip, a.fVars[i].n_rows, a.fVars[i].n_cols );
  	fVars[i] = a.fVars[i];
	}

};

//---------------------------------------------------------------------

void APS_VariableT::operator +=  (const double &a)
{
for (int i=0; i<n_vars; i++) fVars[i] += a;
};


//---------------------------------------------------------------------

void APS_VariableT::operator -=  (const double &a)
{
for (int i=0; i<n_vars; i++) fVars[i] -= a;
};

//---------------------------------------------------------------------

void APS_VariableT::operator *=  (const double &a)
{
for (int i=0; i<n_vars; i++) fVars[i] *= a;
};

//---------------------------------------------------------------------

void APS_VariableT::operator /=  (const double &a) 
{
for (int i=0; i<n_vars; i++) fVars[i] /= a;
};

//---------------------------------------------------------------------

void APS_VariableT::SumOf (APS_VariableT &a,APS_VariableT &b) // Initializes
{
n_vars = a.n_vars;
fVars.Dimension( n_vars );

for (int i=0; i<n_vars; i++) {
  fVars[i].FEA_Dimension( a.fVars[i].n_ip, a.fVars[i].n_rows, a.fVars[i].n_cols );
	fVars[i].SumOf(a.fVars[i],b.fVars[i]); 
}
};

//---------------------------------------------------------------------
