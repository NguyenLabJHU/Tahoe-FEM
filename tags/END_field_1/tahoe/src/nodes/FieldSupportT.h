/* $Id: FieldSupportT.h,v 1.1.2.1 2002-04-30 00:07:12 paklein Exp $ */
#ifndef _FIELD_SUPPORT_T_H_
#define _FIELD_SUPPORT_T_H_

/* forward declarations */
class FEManagerT;
class ElementMatrixT;
template <class TYPE> class nArrayT;
class dArrayT;
class ifstreamT;
class ofstreamT;

/** support for FieldT. Limited interface to get information out
 * of a FieldT. Wrapper for functions in FEManagerT. */
class FieldSupportT
{
public:

	/** constructor */
	FieldSupportT(const FEManagerT& fe);

	/** \name assembly functions */
	/*@{*/
	void AssembleLHS(int group, const ElementMatrixT& elMat, const nArrayT<int>& eqnos) const;
	void AssembleRHS(int group, const dArrayT& elRes, const nArrayT<int>& eqnos) const;
	/*@}*/
	
	/** \name streams */
	/*@{*/
	ifstreamT& Input(void) const;
	ofstreamT& Output(void) const;
	/*@}*/

private:

	/** the boss */
	const FEManagerT& fFEManager;
};

/* constructor */
inline FieldSupportT::FieldSupportT(const FEManagerT& fe):
	fFEManager(fe)
{

}

#endif /* _FIELD_SUPPORT_T_H_ */
