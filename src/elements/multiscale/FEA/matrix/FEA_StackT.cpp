/* $Id: FEA_StackT.cpp,v 1.3 2003-02-03 04:40:24 paklein Exp $ */
#include "FEA.h"

using namespace Tahoe;

//------------------------------------

FEA_StackT::FEA_StackT(void)  							// Constructor  
{ 
	int i,l,stack_length=100,max_ip=100;

	Stack.Dimension(stack_length); 						// Two storage areas (alternate between the two) 
	Shallow_Stack.Dimension(stack_length); 		

	for (i=0; i<stack_length; i++) {
		Stack[i].Allocate(max_ip);           		
		Shallow_Stack[i].Allocate(max_ip);      
		for (l=0; l<max_ip; l++)        				
			Stack[i].vec_ptrs[l] = new double;   	// Allocate space to deep stack to deep copy to
  }

	index  = 0;
	index2 = 0;
}	

//------------------------------------

/* destructor */
FEA_StackT::~FEA_StackT(void)
{
	for (int i = 0; i < Stack.Length(); i++)
		for (int l = 0; l < Stack[i].length; l++)        				
			delete Stack[i].vec_ptrs[l];
}
