/* $Id: DiffusionMatSupportT.h,v 1.1.2.1 2002-10-28 06:49:15 paklein Exp $ */
#ifndef _DIFF_MAT_SUPPORT_T_H_
#define _DIFF_MAT_SUPPORT_T_H_

/* base class */
#include "MaterialSupportT.h"

/* direct members */
#include "dArrayT.h"
#include "dMatrixT.h"

namespace Tahoe {

/* forward declarations */
class LocalArrayT;

/** support for the finite strain Tahoe materials classes */
class DiffusionMatSupportT: public MaterialSupportT
{
public:

	/** constructor */
	DiffusionMatSupportT(int nsd, int ndof, int nip);

	/** destructor */
	~DiffusionMatSupportT(void);

	/** \name field gradients */
	/*@{*/
	/** field gradient at the current integration point */
	const dArrayT& Gradient(void) const;

	/** field gradient at the specified integration point */
	const dArrayT& Gradient(int ip) const;
	/*@}*/

  private:

  	/** \name work space  */
  	/*@{*/
  	ArrayT<dArrayT>  fG_List; /**< field gradient */
  	dArrayT          G_all;   /**< grouped memory for all deformation gradients */
  	/*@}*/	
};

} /* namespace Tahoe */
#endif /* _DIFF_MAT_SUPPORT_T_H_ */
