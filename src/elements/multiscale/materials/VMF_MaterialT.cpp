// $Id: VMF_MaterialT.cpp,v 1.2 2003-02-03 04:40:28 paklein Exp $
#include "VMF_MaterialT.h"

using namespace Tahoe;

VMF_MaterialT::VMF_MaterialT ( void ) { }

int& 		VMF_MaterialT::Number 	( void ) 									{ return matl_no; 					}
void		VMF_MaterialT::Number 	( int i)									{ matl_no = i; 							}
void  	VMF_MaterialT::Assign 	( int index,double val ) 	{ Parameter[index] = val; 	}
void  	VMF_MaterialT::Assign 	( double *a ) 						{ Parameter.Copy (a); 			}
double&	VMF_MaterialT::Retrieve ( int index ) 						{ return Parameter[index]; 	}

//--------------------------------------------------------------

void VMF_MaterialT::Print ( void )
{ 
	for (int i=0; i<n_mp; i++)
		cout << "... INFO >> parameter "<<i<<" = "<<Parameter[i]<<"\n";
	cout << "\n\n";
}
			
//----------------------------------------------------------------

void VMF_MaterialT::E_Nu_2_Lamda_Mu ( int iE,int iNu,int iLamda,int iMu )
{
	double E 	= Parameter[iE];
	double nu = Parameter[iNu];

	if ( E==0 || nu==0) 
		cout << "...WARNING >> VMF_MaterialT::E_Nu_2_Lamda_Mu () : E or nu = 0 \n";	

	Parameter[iLamda] = nu*E/( (1.0+nu)*(1.0-2.0*nu) );
	Parameter[iMu]    = E/( 2*(1.0+nu) );

}

