// $Id: VMS_VariableT.cpp,v 1.3 2003-02-03 04:40:25 paklein Exp $
#include "FEA.h"
#include "VMS.h"

//---------------------------------------------------------------------
/** constructor */

using namespace Tahoe;

VMS_VariableT::VMS_VariableT (const FEA_dMatrixT& GRAD_ua, const FEA_dMatrixT& GRAD_ub) 
{
	Construct	(GRAD_ua, GRAD_ub);
}

//---------------------------------------------------------------------
/** Data initialization: Allocate space for Fa, Fb, GRAD_ua, GRAD_ub  */

void VMS_VariableT::Construct (const FEA_dMatrixT& GRAD_ua, const FEA_dMatrixT& GRAD_ub) 
{
  n_vars = VMS::kNUM_VMS_VARS;  
  fVars.Dimension( n_vars );  

	fVars[VMS::kGRAD_ua] = GRAD_ua; // This = opr allocates if LHS Length=0
	fVars[VMS::kGRAD_ub] = GRAD_ub; 
	fVars[VMS::kGRAD_u].SumOf(GRAD_ua, GRAD_ub);  // This SumOf() checks and allocates if LHS Length=0 
}

//----------------------------------------------------

void VMS_VariableT::Delete_Vars	( void )
{
	for (int i=0; i<n_vars; i++) 
		fVars[i].FEA_Delete(); // ArrayT checks if fLength=0 before deleting
}

//----------------------------------------------------

void VMS_VariableT::Print() { Print(" "); } 

void VMS_VariableT::Print(char *c) { // overload << later

  cout <<"\n VMS_VariableT "<< c <<" follows: \n\n"; 

	for (int l=0; l<n_vars; l++) 
		if (fVars[l].n_ip==0)
			cout << "VMS_VariableT n["<<l<<"] Unallocated \n\n";
    else {
  		cout << "\n Matrix "<<l<<" evaluated at "<<fVars[l].n_ip<<" inegration points (ip): \n"; 
  		for (int i=0; i<fVars[l].n_ip; i++) 
   	 		cout <<"\n "<< c <<" @ ip "<<i<<": \n\n"<< fVars[l][i] << "\n";
			cout << "\n";
		}
	
}

//---------------------------------------------------------------------
//** Retrieve/Fetch/Get either F, F^-1, Fa, Fb, GRAD_ua, or GRAD_ub from class work space
const FEA_dMatrixT& VMS_VariableT::Get(VMS::VarT variable)  
{

  if (!fVars[variable].Length())  // 0 rows indicated un-allocation
    Allocate_and_Compute_Variables(variable);

  return fVars[variable];
  
}

//---------------------------------------------------------------------
/** Compute and store F^-1, Fa, Fb, GRAD(ua), GRAD(ub) : recursive routine */ 
void VMS_VariableT::Allocate_and_Compute_Variables(VMS::VarT kVariable)
{

    switch (kVariable) {

      case VMS::kgrad_u : // grad_u 
        if (!fVars[VMS::kFi].Length())  
          Allocate_and_Compute_Variables(VMS::kFi);     
        fVars[VMS::kgrad_u].MultAB( fVars[VMS::kGRAD_u], fVars[VMS::kFi] ); 
				break;

      case VMS::kgrad_ua : // grad_ua 
        if (!fVars[VMS::kFi].Length())  
          Allocate_and_Compute_Variables(VMS::kFi);     
        fVars[VMS::kgrad_ua].MultAB( fVars[VMS::kGRAD_ua], fVars[VMS::kFi] ); 
				break;

      case VMS::kgrad_ub : // grad_ub 
        if (!fVars[VMS::kFi].Length())  
          Allocate_and_Compute_Variables(VMS::kFi);     
        fVars[VMS::kgrad_ub].MultAB( fVars[VMS::kGRAD_ub], fVars[VMS::kFi] ); 
				break;

      case VMS::kF :  // F
        fVars[VMS::kF] = fVars[VMS::kGRAD_u];
				fVars[VMS::kF].PlusIdentity();
				break;

      case VMS::kFi : // F^-1 
        if (!fVars[VMS::kF].Length())  
          Allocate_and_Compute_Variables(VMS::kF);     
        fVars[VMS::kFi] = fVars[VMS::kF];
        fVars[VMS::kFi].Inverse();                            
				break;
    
			case VMS::kFa : // Fa := F^Alpha
        if (!fVars[VMS::kGRAD_ua].Length())  
          Allocate_and_Compute_Variables(VMS::kGRAD_ua);     
        fVars[VMS::kFa] = fVars[VMS::kGRAD_ua];
        fVars[VMS::kFa].PlusIdentity();      
				break;
    
      case VMS::kFai : // (F^a)^-1
        if (!fVars[VMS::kFa].Length())  
          Allocate_and_Compute_Variables(VMS::kFa);     
        fVars[VMS::kFai] = fVars[VMS::kFa];
        fVars[VMS::kFai].Inverse();         
				break;
    
			case VMS::kFb : // Fb := F^Beta = 1 +  [d(ub)/dX]Fa^-1
        if (!fVars[VMS::kGRAD_ub].Length())  
          Allocate_and_Compute_Variables(VMS::kGRAD_ub);     
        if (!fVars[VMS::kFai].Length())  
          Allocate_and_Compute_Variables(VMS::kFai);     
        fVars[VMS::kFb].MultAB( fVars[VMS::kGRAD_ub], fVars[VMS::kFai] );
        fVars[VMS::kFb].PlusIdentity();                     
				break;

      case VMS::kFbi : // Fb^-1
        if (!fVars[VMS::kFb].Length())  
          Allocate_and_Compute_Variables(VMS::kFb);     
        fVars[VMS::kFbi] = fVars[VMS::kFb];
        fVars[VMS::kFbi].Inverse();                    
				break;

      default:
        cout << "\n VMS_VariableT:::Allocate_and_Compute_Variables() bad VarT type" << endl;
				//throw eGeneralFail;
				break;
    }
}

//---------------------------------------------------------------------

void VMS_VariableT::operator  =  (const VMS_VariableT &a)	// Initializes
{
	n_vars = a.n_vars;
	fVars.Dimension( n_vars );

	for (int i=0; i<n_vars; i++) {
  	fVars[i].FEA_Dimension( a.fVars[i].n_ip, a.fVars[i].n_rows, a.fVars[i].n_cols );
  	fVars[i] = a.fVars[i];
	}

};

//---------------------------------------------------------------------

void VMS_VariableT::operator +=  (const double &a)
{
for (int i=0; i<n_vars; i++) fVars[i] += a;
};


//---------------------------------------------------------------------

void VMS_VariableT::operator -=  (const double &a)
{
for (int i=0; i<n_vars; i++) fVars[i] -= a;
};

//---------------------------------------------------------------------

void VMS_VariableT::operator *=  (const double &a)
{
for (int i=0; i<n_vars; i++) fVars[i] *= a;
};

//---------------------------------------------------------------------

void VMS_VariableT::operator /=  (const double &a) 
{
for (int i=0; i<n_vars; i++) fVars[i] /= a;
};

//---------------------------------------------------------------------

void VMS_VariableT::SumOf (VMS_VariableT &a,VMS_VariableT &b) // Initializes
{
n_vars = a.n_vars;
fVars.Dimension( n_vars );

for (int i=0; i<n_vars; i++) {
  fVars[i].FEA_Dimension( a.fVars[i].n_ip, a.fVars[i].n_rows, a.fVars[i].n_cols );
	fVars[i].SumOf(a.fVars[i],b.fVars[i]); 
}
};

//---------------------------------------------------------------------
