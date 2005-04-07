/* $Id: ABAQUS_BaseT.h,v 1.3 2004-01-05 23:36:01 paklein Exp $ */
#ifndef _ABAQUS_BASE_T_H_
#define _ABAQUS_BASE_T_H_

/* library support options */
#ifdef __F2C__

/* f2c */
#include "ArrayT.h"
#include "f2c.h"
//MIPSpro 7.2.1 chokes on f2c.h unless some C++ headers are read first
//compaq-cxx-6.5.9.28 also chokes on f2c.h and needs a more extensive
//C++ header to precede f2c.h

namespace Tahoe {

/* forward declarations */
class ifstreamT;
template <class nTYPE> class nMatrixT;
class dMatrixT;
class dSymMatrixT;
class ElementCardT;
class StringT;
template <class nTYPE> class nArrayT;

/** some ABAQUS basics */
class ABAQUS_BaseT
{
public:

	/** constructor */
	ABAQUS_BaseT(void);

protected:

	/** name stress conversion functions */
	/*@{*/
	void dMatrixT_to_ABAQUS(const dMatrixT& A, nMatrixT<doublereal>& B) const;
	void ABAQUS_to_dSymMatrixT(const doublereal* pA, dSymMatrixT& B, bool convert_shear = false) const;
	void dSymMatrixT_to_ABAQUS(const dSymMatrixT& A, doublereal* pB, bool convert_shear = false) const;
	/*@}*/

	/** \name read ABAQUS-format input */
	/*@{*/
	void Read_ABAQUS_Input(ifstreamT& in, StringT& name, nArrayT<doublereal>& properties,
		double& density, integer& nstatv, bool& nonsym) const;
	bool Next_ABAQUS_Keyword(ifstreamT& in) const;
	bool Skip_ABAQUS_Symbol(ifstreamT& in, char c) const;
	void Skip_ABAQUS_Comments(ifstreamT& in) const;
	void Read_ABAQUS_Word(ifstreamT& in, StringT& word, bool to_upper = true) const;
	/*@}*/	
};

} /* namespace Tahoe */

#else /* __F2C__ */

#ifndef __MWERKS__
#error "ABAQUS_BaseT requires __F2C__"
#endif

#endif /* __F2C__ */

#endif /* _ABAQUS_BASE_T_H_ */
