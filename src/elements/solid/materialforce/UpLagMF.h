/* $Id: UpLagMF.h,v 1.4 2003-11-12 19:21:19 thao Exp $ */

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
    virtual void FormKd(double constK);

 protected:
    /*material force evaluation*/
    void ComputeMatForce(dArray2DT& output);
    void MatForceVolMech(dArrayT& elem_val);
    void MatForceDissip(dArrayT& elem_val, const dArray2DT& internalstretch);
    void MatForceSurfMech(dArrayT& global_val);

    /*utility funtions*/
    /*extrapolate element ip state variables values to nodes*/
    void Extrapolate(void);
  
 protected:	

    /*current material*/
    FSSolidMatT* fCurrFSMat;

    /*connectivities*/
    RaggedArray2DT<int> fXConnects;

    /*localization check*/
    iArrayT floc_flags;
    dArray2DT felem_centers;
    dArray2DT fnormals;
    
 private:
    dMatrixT fEshelby;        /*eshelby energy momentum tensor*/
    dSymMatrixT fC;
    /*nodal and interpolated body force*/
    dArrayT fBodyForce;       /*body and inertial force vector*/
    dArrayT fip_body;         /*body force at integration point*/
    
    /*material force*/
    dArrayT fMatForce;        /*nodal material force vector*/
    dArrayT fDissipForce;     /*nodal dissipation force vector*/
    dArrayT felem_rhs;        /*element force vector*/

    /*internal variables*/
    iArrayT fInternalDOF;     /*dof of internal variable tensors*/
    int fNumInternalVal;      /*total number of internal variables*/

    dArrayT fGlobalMass;      /*lumped global mass matrix*/
    dArray2DT fGlobalVal;     /*global extrapolation values*/

    dArrayT felem_mass;       /*lumped element mass matrix*/
    dArray2DT felem_val;      /*element extrapolation values*/
    
    dArray2DT fGradInternalStrain;  /*gradient of internal strains at ip*/

    /*surface variables*/
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
