/* $Id: SSMF.h,v 1.3 2003-11-19 06:09:46 thao Exp $ */

#ifndef _SSMF_H_
#define _SSMF_H_

/* base class */
#include "SmallStrainT.h"
#include "MFSupportT.h"
#include "ofstreamT.h"
#include "RaggedArray2DT.h"
namespace Tahoe {

/* forward declarations */
class SSSolidMatT;
class SSMatSupportT;
class ifstreamT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class SSMF: public SmallStrainT, public MFSupportT
{
  public:
    /** constructor */
    SSMF(const ElementSupportT& support, const FieldT& field);

    virtual void Initialize(void);
    virtual void SetGlobalShape(void);
    
    /*accessor for displacement gradient*/
    const dMatrixT& DisplacementGradient(void) const;

    /*output*/
    /* register self for output */
    virtual void RegisterOutput(void);

    /* send output */
    virtual void WriteOutput(void);

    /*connectivities for parrallel computing*/
    virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
			   AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
 
 protected:
    /*material force evaluation*/
    void ComputeMatForce(dArray2DT& output);
    void MatForceVolMech(dArrayT& elem_val);
    void MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch);
    void MatForceDynamic(dArrayT& elem_val);
    void MatForceSurfMech(dArrayT& global_val);

    /*utility funtions*/
    /*extrapolate element ip state variables values to nodes*/
    void Extrapolate(void);
  
 protected:	

    /*current material*/
    SSSolidMatT* fCurrSSMat;

    /*connectivities*/
    RaggedArray2DT<int> fXConnects;

 private:
    ArrayT<dMatrixT> fGradU_List;
    dMatrixT fEshelby;

    /*nodal and interpolated body force*/
    dArrayT fBodyForce;
    dArrayT fip_body;
    
    /*material force */
    dArrayT fMatForce;
    dArrayT fDissipForce;
    dArrayT fDynForce;
    dArrayT felem_rhs;

    /*dynamic analysis variables*/
    bool fdynamic;           /*flag for dynamic analysis*/
    dArrayT fVel;             /*integration point velocity vector*/
    dArrayT fAcc;             /*integration point acceleration vector*/
    dMatrixT fGradVel;        /*integration point velocity gradient*/

    /*internal variables*/
    iArrayT fInternalDOF;
    int fNumInternalVal;

    dArrayT fGlobalMass;
    dArray2DT fGlobalVal;

    dArrayT felem_mass;
    dArray2DT felem_val;
    
    dArray2DT fGradInternalStrain;

    /*crack surface evaluation*/
    LocalArrayT ftraction;    /*traction at element facet*/
    LocalArrayT fsurf_disp;   /*displacement at element facet*/
    LocalArrayT fsurf_coords; /*coordinates of element facet*/

    dArrayT fsurf_val;        /*element force vector at facet*/
    dArrayT fip_tract;        /*traction at element facet integration points*/
    dMatrixT fgradU;          /*displacement gradient along facet direction at facet integration points*/
    dMatrixT fjacobian;       /*surface jacobian mapping global to parent coordinate*/
    dMatrixT fjac_loc;        /*surface jacobian mapping local to parent coordinates*/
    dMatrixT fQ;              /*rotation from local to global coords*/
};

/* inlines */
inline const dMatrixT& SSMF::DisplacementGradient(void) const
{
	return fGradU_List[CurrIP()];
}

} // namespace Tahoe 
#endif /* _SSMF_H_ */
