/* $Id: FiniteStrainMF.h,v 1.1 2003-02-12 18:37:41 thao Exp $ */
#ifndef _FiniteStrain_MF_H_
#define _FiniteStrain_MF_H_

/* base class */
#include "FiniteStrainT.h"

namespace Tahoe {

/* forward declarations */
class FSMatSupportT;
class FSSolidMatT;

/** Interface for linear strain deformation and field gradients */
class FiniteStrainMF: public FiniteStrainT
{
 public:
      
	/* constructor */
	FiniteStrainMF(const ElementSupportT& support, const FieldT& field);
	/*destructor*/
        ~FiniteStrainMF(void);

	/* register self for output */
	virtual void RegisterOutput(void);
	/* send output */
	virtual void WriteOutput(void);

 protected:
	/*material force evaluation*/
	void ComputeMatForce(dArray2DT& output);
     
	void MatForceVolMech(dArrayT& elem_val);
	void MatForceDissip(dArrayT& elem_val, const dArrayT& statev);
	void MatForceSurfMech(dArrayT& global_val);

	/*Assemble nodal material force vectors for element group*/
	void AssembleMatForce(const dArrayT& elem_val, dArrayT& global_val,
			      const iArrayT& nodes);
	/* map nodal ordering of element group*/ 
	void MapOutput(void);

 protected:
	/*current material*/
	FSSolidMatT* fCurrFSMat;
	
	/* material force output ID */
	int fMatForceOutputID;
	OutputSetT fOutputSet;
	
	/*output dimension*/
	int fNumGroupNodes;
	dArrayT fMap;

};

} // namespace Tahoe 
#endif /* _FINITE_STRAIN_MF_H_ */
