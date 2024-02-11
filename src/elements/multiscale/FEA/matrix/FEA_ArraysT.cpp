#include "FEA.h"

using namespace Tahoe;

//############################################################
//######### MATRIX ARRAY #####################################
//############################################################

FEA_dMatrix_ArrayT::FEA_dMatrix_ArrayT (int n_mat,int n_ip,int n_rows,int n_cols) : ArrayT <FEA_dMatrixT> (n_mat) 
{
	for (int i=0; i<fLength; i++) 	(*this)[i].FEA_Dimension (n_ip, n_rows, n_cols); 
}	

//----------------------------------------------------

void FEA_dMatrix_ArrayT::Construct (int n_mat,int n_ip,int n_rows,int n_cols) 
{
	Dimension ( n_mat );
	for (int i=0; i<fLength; i++) 	(*this)[i].FEA_Dimension (n_ip, n_rows, n_cols); 
}	

//----------------------------------------------------

void FEA_dMatrix_ArrayT::Print() { Print(" "); } 

void FEA_dMatrix_ArrayT::Print(char *c) { // overload << later

  char string[10];

  cout <<"\n FEA_dMatrix_ArrayT "<< c <<" follows: \n\n"; 

	for (int l=0; l<fLength; l++) 
		if ((*this)[l].n_ip==0)
			cout << " FEA_dScalar_ArrayT: n["<<l<<"] Unallocated \n\n";
    else {
			sprintf(string,"%d",l);
			(*this)[l].Print(string);
		}
}

//############################################################
//######### VECTOR ARRAY #####################################
//############################################################

FEA_dVector_ArrayT::FEA_dVector_ArrayT (int n_vec,int n_ip,int n_rows) : ArrayT <FEA_dVectorT> (n_vec) 
{
	for (int i=0; i<fLength; i++) 	(*this)[i].FEA_Dimension (n_ip, n_rows); 
}	

//----------------------------------------------------

void FEA_dVector_ArrayT::Construct (int n_vec,int n_ip,int n_rows) 
{
	Dimension ( n_vec );
	for (int i=0; i<fLength; i++) 	(*this)[i].FEA_Dimension (n_ip, n_rows); 
}	

//----------------------------------------------------

void FEA_dVector_ArrayT::Print() { Print(" "); } 

void FEA_dVector_ArrayT::Print(char *c) { // overload << later

  char string[10];

  cout <<"\n FEA_dVector_ArrayT "<< c <<" follows: \n\n"; 

	for (int l=0; l<fLength; l++) 
		if ((*this)[l].Length()==0)
			cout << " FEA_dScalar_ArrayT: n["<<l<<"] Unallocated \n\n";
    else {
			sprintf(string,"%d",l);
			(*this)[l].Print(string);
		}
}

//############################################################
//######### SCALAR ARRAY #####################################
//############################################################

FEA_dScalar_ArrayT::FEA_dScalar_ArrayT (int n_scal,int n_ip) : ArrayT <FEA_dScalarT> (n_scal) 
{
	for (int i=0; i<fLength; i++) 	(*this)[i].Dimension (n_ip); 
}	

//----------------------------------------------------

void FEA_dScalar_ArrayT::Construct (int n_scal,int n_ip) 
{
	Dimension ( n_scal );
	for (int i=0; i<fLength; i++) 	(*this)[i].Dimension (n_ip); 
}	

//----------------------------------------------------

void FEA_dScalar_ArrayT::Print() { Print(" "); } 

void FEA_dScalar_ArrayT::Print(char *c) { // overload << later

  char string[10];

  cout <<"\n FEA_dScalar_ArrayT "<< c <<" follows: \n\n"; 

	for (int l=0; l<fLength; l++) 
		if ((*this)[l].Length()==0)
			cout << " FEA_dScalar_ArrayT: n["<<l<<"] Unallocated \n\n";
    else {
			sprintf(string,"%d",l);
			(*this)[l].Print(string);
		}
}

