/* $Id: UpLagMF.h,v 1.7 2003-11-24 17:35:13 thao Exp $ */

#ifndef _UpLagMF_H_
#define _UpLagMF_H_

/* base class */
#include "UpdatedLagrangianT.h"
#include "MFSupportT.h"
#include "LocalizeT.h"
#include "ofstreamT.h"
#include "RaggedArray2DT.h"
namespace Tahoe {

/* forward declarations */
class FSSolidMatT;
class ifstreamT;
class StringT;

/** Interface for linear strain deformation and field gradients */
class UpLagMF: public UpdatedLagrangianT, public MFSupportT, public LocalizeT
{
  public:
    /** constructor */
    UpLagMF(const ElementSupportT& support, const FieldT& field);

    virtual void Initialize(void);
    virtual void SetGlobalShape(void);

    /*output*/
    /* register self for output */
    virtual void RegisterOutput(void);

    /* send output */
    virtual void WriteOutput(void);

    /*connectivities for parrallel computing*/
    virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
			   AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;

    /*form right hand side*/
    virtual GlobalT::RelaxCodeT RelaxSystem(void);

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
    FSSolidMatT* fCurrFSMat;

    /*connectivities*/
    RaggedArray2DT<int> fXConnects;

 private:
    dMatrixT fEshelby;        /*eshelby energy momentum tensor*/
    dSymMatrixT fC;

    /*material force*/
    dArrayT fMatForce;        /*nodal material force vector*/
    dArrayT fDissipForce;     /*nodal dissipation force vector*/
    dArrayT fDynForce;
    dArrayT felem_rhs;        /*element force vector*/

    /*nodal and interpolated body force*/
    dArrayT fBodyForce;       /*nodal body force vector*/
    dArrayT fip_body;         /*body force at integration point*/
 
    /*dynamic analysis variables*/
    bool fdynamic;           /*flag for dynamic analysis*/
    dArrayT fVel;             /*integration point velocity vector*/
    dArrayT fAcc;             /*integration point acceleration vector*/
    dMatrixT fGradVel;        /*integration point velocity gradient*/

    /*internal variables for inelastic materials*/
    iArrayT fInternalDOF;     /*dof of internal variable tensors*/
    int fNumInternalVal;      /*total number of internal variables*/

    dArrayT fGlobalMass;      /*lumped global mass matrix*/
    dArray2DT fGlobalVal;     /*global extrapolation values*/

    dArrayT felem_mass;       /*lumped element mass matrix*/
    dArray2DT felem_val;      /*element extrapolation values*/
    
    dArray2DT fGradInternalStrain;  /*gradient of internal strains at ip*/

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

} // namespace Tahoe 
#endif /* _UpLagMF_H_ */
