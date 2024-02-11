// $Id: FEA_StackT.h,v 1.4 2003-02-03 04:40:24 paklein Exp $
#ifndef _FEA_STACKT_H_
#define _FEA_STACKT_H_

namespace Tahoe {

class FEA_StackT 
{
  friend class FEA_dMatrixT;             // a  *  b   *  c
  friend class FEA_dVectorT;             //  \   / \    /
  friend class FEA_dScalarT;             //   \ /   \  /
  friend class FEA_EquateT;              //    d  *   e
                                         //     \    /
    public:                              //      \  /
                                         //       f    
  	FEA_StackT(void);

	/** destructor */
	~FEA_StackT(void);

		/** Circular, alternates 0,1,0,1,0,1,... */

 	 	//bool Next_Stack(void) 				{ index  = (index==0)  ? 1 : 0; return (index); 	} 
 	 	//bool Next_Shallow_Stack(void) { index2 = (index2==0) ? 1 : 0; return (index2); 	} 

 	 	int Next_Stack (void) 				{ return index++;  					}
 	 	int Next_Shallow_Stack (void) { return index2++; 					} 
	  void Reset (void) 						{	index = index2 = 0; }

	private:
   
  	ArrayT<FEA_EquateT> Stack; 					// Stack with deep copies (new values due to arithmatic results written here)
  	ArrayT<FEA_EquateT> Shallow_Stack; 	// Stack with shallow copies (copies of pointers to components copied here)

		int index;
		int index2;

		//bool index;
		//bool index2;

		// ex.  A(1,2) = B(3,4)*C(5,6);
		// NOTE: LHS needs to copy shallow pointers so data will change, RHS components will be shallow pointers,
		// but there product will be deep (new data)
};

}
#endif

