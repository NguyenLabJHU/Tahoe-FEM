// $Id: APS_MaterialT.cpp,v 1.2 2004-02-04 00:40:52 raregue Exp $

#include "APS_MaterialT.h"

using namespace Tahoe;

APS_MaterialT::APS_MaterialT ( void ) { }

int& 		APS_MaterialT::Number 	( void ) 					{ return matl_no; 			}
void		APS_MaterialT::Number 	( int i)					{ matl_no = i; 				}
void  		APS_MaterialT::Assign 	( int index,double val ) 	{ Parameter[index] = val; 	}
void  		APS_MaterialT::Assign 	( double *a ) 				{ Parameter.Copy (a); 		}
double&		APS_MaterialT::Retrieve ( int index ) 				{ return Parameter[index]; 	}

//--------------------------------------------------------------

void APS_MaterialT::Print ( void )
{ 
	for (int i=0; i<n_mp; i++)
		cout << "... INFO >> parameter "<<i<<" = "<<Parameter[i]<<"\n";
	cout << "\n\n";
}


