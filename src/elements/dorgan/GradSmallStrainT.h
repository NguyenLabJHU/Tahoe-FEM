/* $Id: GradSmallStrainT.h,v 1.1 2003-11-19 20:38:22 rdorgan Exp $ */ 
#ifndef _GRAD_SMALL_STRAIN_T_H_ 
#define _GRAD_SMALL_STRAIN_T_H_ 

/* base classes */
#include "SmallStrainT.h"

/* direct members */
#include "LocalArrayT.h"
//#include "GeometryT.h"
#include "C1GeometryT.h"
#include "dSymMatrixT.h"

namespace Tahoe {

/* forward declarations */
class GradSSSolidMatT;
class GradSSMatSupportT;
//class SolidMatListT;
class C1ShapeFunctionT;

/** element formulation with gradient plasticity constitutive model */
class GradSmallStrainT: public SmallStrainT
{
public:

        /** constructor */
        GradSmallStrainT(const ElementSupportT& support, const FieldT& disp, 
                             const FieldT& iso_hard);

        /** destructor */
        ~GradSmallStrainT(void);
        
        /** data initialization */
        virtual void Initialize(void); 
        
        /** collecting element group equation numbers. See ElementBaseT::Equations
         * for more information */
        virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
                               AutoArrayT<const RaggedArray2DT<int>*>& eq_2);
        
        /** \name isotropic hardening */
        /*@{*/
        double& LinearR(void) const { return fR_List[CurrIP()]; };
        double& LinearR(int ip) const { return fR_List[ip]; };
        /*@}*/
        
        /** \name isotropic hardening from the end of the previous time step */
        /*@{*/
        double& LinearR_last(void) const { return fR_last_List[CurrIP()]; };
        double& LinearR_last(int ip) const { return fR_last_List[ip]; };
        /*@}*/
        
        /** \name Laplacian isotropic hardening */
        /*@{*/
        double& LinearLaplacianR(void) const { return fLapR_List[CurrIP()]; };
        double& LinearLaplacianR(int ip) const { return fLapR_List[ip]; };
        /*@}*/
        
        /** \name Laplacian isotropic hardening from the end of the previous time step */
        /*@{*/
        double& LinearLaplacianR_last(void) const { return fLapR_last_List[CurrIP()]; };
        double& LinearLaplacianR_last(int ip) const { return fLapR_last_List[ip]; };
        /*@}*/
        
        /** return the number of degrees of freedom for iso_hard per node */
        int NumDOF_R(void) const { return fIsoHardening.NumDOF();} ;
        
        /** number of element integration points for iso_hard field */
        int NumIP_R(void) const { return fNumIP_R;} ;
        
        /** reference to element shape functions */
        const C1ShapeFunctionT& C1ShapeFunction(void) const;
        
        /** return the geometry code */
        C1GeometryT::CodeT GeometryCode_R(void) const;
        
protected:
        
        /** construct a new material support and return a pointer. Recipient is responsible for
         * for freeing the pointer.
         * \param p an existing MaterialSupportT to be initialized. If NULL, allocate
         *        a new MaterialSupportT and initialize it. */
        virtual MaterialSupportT* NewMaterialSupport(MaterialSupportT* p = NULL) const;
        
        /** initialization functions */
        virtual void SetLocalArrays(void);
        virtual void SetShape(void);

        /** form shape functions and derivatives */
        virtual void SetGlobalShape(void);
        
        /** form the element stiffness matrix */
        virtual void FormStiffness(double constK);
        
        /** calculate the internal force contribution ("-k*d") */
        virtual void FormKd(double constK);
        
        /** increment current element */
        virtual bool NextElement(void);        
        
        /** return a const reference to the run state flag */
        virtual GlobalT::SystemTypeT TangentType(void) const;
        
        /** write element group parameters to out */
        virtual void PrintControlData(ostream& out) const;
        
private:
        /** set the \e h matrix using the given shape functions */
        virtual void Set_h(dMatrixT& h) const;
        
        /** set the \e p matrix using the given shape functions */
        virtual void Set_p(dMatrixT& p) const;

protected:
        /** \name return values */
        /*@{*/
        dArrayT fR_List;
        dArrayT fR_last_List;
        
        dArrayT fLapR_List;
        dArrayT fLapR_last_List;
        /*@}*/
          
        /** \name element iso_hard in local ordering for current element */
        /*@{*/
        LocalArrayT fLocR;      /**< hardness */
        LocalArrayT fLocLastR;  /**< hardness from last time increment */

        dArrayT fLocRTranspose;      /**< hardness */
        dArrayT fLocLastRTranspose;  /**< hardness from last time increment */
        /*@}*/
        
private:
        /* \name fields */
        /*@{*/
        const FieldT& fDisplacement; /**< displacement field */
        const FieldT& fIsoHardening; /**< hardening parameter fiel */
        /*@}*/
        
        /** \name shape functions for isotropic hardening */
        C1ShapeFunctionT* fShapes_R;
        
        /** the material support used to construct materials lists. This pointer
         * is only set the first time GradSmallStrainT::NewMaterialList is called. */
        GradSSMatSupportT* fGradSSMatSupport;
        
        /* run time */
        GradSSSolidMatT*  fCurrMaterial_Grad;
        
        /** \name work space */
        /*@{*/
        /** shape functions for R */
        dMatrixT fh;      /**< C1 shape functions */
        dMatrixT fhT;     /**< C1 shape functions (Transpose) */
        
        /** shape functions for LapR */
        dMatrixT fp;      /**< Laplacian of C1 shape functions */
        dMatrixT fpT;     /**< Laplacian of C1 shape functions (Transpose) */
        
        /** stiffnesses */
        ElementMatrixT fKaa;         /**< elastic stiffness matrix */
        ElementMatrixT fKar, fKra;   /**< off-diagonal matrices */
        ElementMatrixT fKrr_hh;      /**< local part of Krr (Krr: gradient dependent diagonal matrix) */
        ElementMatrixT fKrr_hh_constraint;      /**< constraint matrix */
        ElementMatrixT fKrr_hp;      /**< gradient part of Krr (Krr: gradient dependent diagonal matrix) */
       
        /** returned matrices */
        dMatrixT fODM_ar;
        dMatrixT fODM_ra;
        dMatrixT fGM;
        dMatrixT fI;
        /*@}*/
        
        /** \name dimensions */
        /*@{*/
        int fNumIP;             /**< number of integration points */
        int fNumIP_R;           /**< number of integration points for iso_hard field */
        
        int fNumSD;             /**< number of spatial dimensions */
        
        //        int fNumElementNodes;   /**< number of element nodes */
        int fNumElements;       /**< number of integration points */
        
        int fNumDOF;            /**< number of degrees of freedom for displacement field */
        int fNumDOF_R;          /**< number of degrees of freedom for iso_hard field */
        int fNumDOF_Total;      /**< number of total degrees of freedom (ndf_Disp + ndf_R) */
        /*@}*/
        
        /* print debug information */
        bool print_GlobalShape;
        bool print_Kd;
        bool print_Stiffness;
        
        /** constraint constant */
        double fKConstraintA;
        double fKConstraintB;
        double fRConstraintA;
        double fRConstraintB;
        
        /** element parameter */
        C1GeometryT::CodeT fGeometryCode_R;
};

/* inlines */

/* return the geometry code */
inline C1GeometryT::CodeT GradSmallStrainT::GeometryCode_R(void) const
{ return fGeometryCode_R; }

/* accessors */
inline const C1ShapeFunctionT& GradSmallStrainT::C1ShapeFunction(void) const
{
#if __option(extended_errorcheck)
        if (!fShapes_R)
        {
                cout << "\n GradSmallStrainT::C1ShapeFunction: no shape functions" << endl;
                throw ExceptionT::kGeneralFail;
        }
#endif
        return *fShapes_R;
}

} // namespace Tahoe 

#endif /* _GRAD_SMALL_STRAIN_T_H_ */
