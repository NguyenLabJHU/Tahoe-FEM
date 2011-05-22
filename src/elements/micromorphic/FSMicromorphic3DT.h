/* $Id: FSMicromorphic3DT.h,v 1.192 2011-05-22 02:32:44 isbuga Exp $ */
//DEVELOPMENT
#ifndef _FS_MICROMORPHIC_3D_T_H_
#define _FS_MICROMORPHIC_3D_T_H_

/* base classes */
#include "ElementBaseT.h"
#include "StringT.h"
#include "Traction_CardT.h"
#include "ShapeFunctionT.h"
#include "eIntegratorT.h"
#include "ModelManagerT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"

#include "iAutoArrayT.h"
#include "ScheduleT.h"

/* direct members */
#include "LocalArrayT.h"
#include "GeometryT.h"
#include "dTensor3DT.h"// I added this header file to use Tensor built-in functions!

#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"


namespace Tahoe {

/* forward declarations */
class ShapeFunctionT;
class Traction_CardT;
class StringT;

/** FSMicromorphic3DT: This class contains a coupled finite deformation
 * micromorphic (displacement and micro-displacement gradient dofs) finite
 * element implementation in 3D.  The model description can be
 * found in Regueiro ASCE JEM 135:178-191 for pressure-sensitive elastoplasticity.
 **/

class FSMicromorphic3DT: public ElementBaseT
{

public:
/* material parameters */
    enum fMaterial_T         {
        kMu,
        kLambda,
        kKappa,
        kNu,
        kSigma_const,
        kTau,
        kEta,
        kRho_0,
        kTau1,
        kTau2,
        kTau3,
        kTau4,
        kTau5,
        kTau6,
        kTau7,
        kTau8,
        kTau9,
        kTau10,
        kTau11,
        // plasticity parameters
        kHc,
        kc0,
        kHc_chi,
        kc0_chi,        
        //kZ0c,
        kFphi,
        kDpsi,
        kFphi_chi,
        kDpsi_chi,
        //
        kg,
        kg1,
        kg2,
        kg3,
        //add to this list
        kNUM_FMATERIAL_TERMS        };

    enum fMaterialState_T         {
    //    kkappa,
        kc,
        kc_chi,            
        //    kZkappa,
    //    kZc,
    //    khkappa,
        khc,
        khc_chi,
    //    kIntrinsic_Perm,
    //    kJ,
    //    kJp,
    //    kphi_s,
    //    kphi_f,
    //    kDevSS,
    //    kMeanS,
    //    kEpsVolp,
        kDelgamma,
        kDelgammachi,
        ktrSigma,
        kSigma_inv,
        ktrRel,
        kRel_inv,
        ktrm,
        km_inv,
        ktrS,
        kinvdevS,
        ktrSIGMA_S,
        kinvdevSIGMA_S,
        kinvtrM,
        kinvdevM,
        kinvPhi,
        kinvGPhi,
        kNUM_FMATERIAL_STATE_TERMS
    };


//  enum fIntegrate_T         {
//      kBeta,
//      kGamma,
//      kNUM_FINTEGRATE_TERMS        };

    /** constructor */
    FSMicromorphic3DT(        const ElementSupportT& support );

    /** destructor */
    ~FSMicromorphic3DT(void);

    /** reference to element shape functions */
    const ShapeFunctionT& ShapeFunctionDispl(void) const;
    const ShapeFunctionT& ShapeFunctionMicro(void) const;

    /** echo input */
    void Echo_Input_Data (void);

    /** return true if the element contributes to the solution of the
     * given group. ElementBaseT::InGroup returns true if group is the
     * same as the group of the FieldT passed in to ElementBaseT::ElementBaseT. */
    virtual bool InGroup(int group) const;

    /* initialize/finalize time increment */
    /*@{*/
    virtual void InitStep(void);
    virtual void CloseStep(void);
    //virtual GlobalT::RelaxCodeT ResetStep(void); // restore last converged state

    /** element level reconfiguration for the current time increment */
    //virtual GlobalT::RelaxCodeT RelaxSystem(void);
    /*@}*/

    /** collecting element group equation numbers. See ElementBaseT::Equations
     * for more information */
    virtual void Equations( AutoArrayT<const iArray2DT*>& eq_d,
                AutoArrayT<const RaggedArray2DT<int>*>& eq_phi);

    /** return a const reference to the run state flag */
    virtual GlobalT::SystemTypeT TangentType(void) const;

    /** accumulate the residual force on the specified node
     * \param node test node
     * \param force array into which to assemble to the residual force */
    virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);

    /** returns the energy as defined by the derived class types */
    virtual double InternalEnergy(void);

    /** \name writing output */
    /*@{*/
    /** register element for output */
    virtual void RegisterOutput(void);

    /** write element output */
    virtual void WriteOutput(void);

    /** compute specified output parameter and send for smoothing */
    virtual void SendOutput(int kincode);
    /*@}*/

    /** return geometry and number of nodes on each facet */
    void FacetGeometry(ArrayT<GeometryT::CodeT>& facet_geometry, iArrayT& num_facet_nodes) const;

    /** return the geometry code */
    virtual GeometryT::CodeT GeometryCode(void) const;

    /*set active elements*/
    //virtual void SetStatus(const ArrayT<ElementCardT::StatusT>& status);

    /** initial condition/restart functions (per time sequence) */
    virtual void InitialCondition(void);

    /** mass types */
    enum MassTypeT {kNoMass = 0, /**< do not compute mass matrix */
            kConsistentMass = 1, /**< variationally consistent mass matrix */
                kLumpedMass = 2, /**< diagonally lumped mass */
             kAutomaticMass = 3  /**< select the mass type base on the time integration scheme */};
    MassTypeT static int2MassTypeT(int i);

    /** \name restart functions */
    /*@{*/
    /** write restart data to the output stream. Should be paired with
     * the corresponding ElementBaseT::ReadRestart implementation. */
    virtual void WriteRestart(ostream& out) const;

    /** read restart data to the output stream. Should be paired with
     * the corresponding ElementBaseT::WriteRestart implementation. */
    virtual void ReadRestart(istream& in);
    /*@}*/

    /** \name implementation of the ParameterInterfaceT interface */
    /*@{*/
    /** information about subordinate parameter lists */
    virtual void DefineSubs(SubListT& sub_list) const;

    /** return the description of the given inline subordinate parameter list */
    virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order,
                     SubListT& sub_lists) const;

    /** a pointer to the ParameterInterfaceT of the given subordinate */
    virtual ParameterInterfaceT* NewSub(const StringT& name) const;

    /** \name implementation of the ParameterInterfaceT interface */
    /*@{*/
    /** describe the parameters needed by the interface */
    virtual void DefineParameters(ParameterListT& list) const;

    /** accept parameter list */
    virtual void TakeParameterList(const ParameterListT& list);
    /*@}*/

protected:

    /** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
    /*@{*/
    /** form group contribution to the stiffness matrix */
    virtual void LHSDriver(GlobalT::SystemTypeT);

    /** form group contribution to the residual */
    virtual void RHSDriver(void);
    /*@}*/

    /** compute shape functions and derivatives */
    virtual void SetGlobalShape(void);

    void Select_Equations ( const int &iBalLinMom, const int &iBalFirstMomMom );

private:

    /** \name solution methods.
     * Both of these drivers assemble the LHS as well as the residual.
     */
    /*@{*/
    /** driver for staggered solution */
    void RHSDriver_staggered(void);

    /** driver for monolithic solution */
    void RHSDriver_monolithic(void);
    /*@}*/

protected:

    /* output control */
    iArrayT        fNodalOutputCodes;
    iArrayT        fElementOutputCodes;

private:

    /** Gradients and other matrices */
    dMatrixT fgrad_u, fgrad_u_n;
    dMatrixT fgrad_Phi, fgrad_Phi_n;

    dMatrixT fShapeDispl, fShapeDisplGrad, fShapeDisplGrad_t;
    dMatrixT fShapeDisplGrad_t_Transpose, fShapeDisplGradGrad, fShapeDisplGrad_temp;
    dMatrixT fShapeMicro,fShapeMicroGrad,fShapeMicroGrad_temp;

    dMatrixT fDefGrad, fDefGradInv, fDefGradInvMatrix;

    /** \name  values read from input in the constructor */
    /*@{*/
    /** element geometry */
    GeometryT::CodeT fGeometryCode_displ, fGeometryCodeSurf_displ;
    GeometryT::CodeT fGeometryCode_micro, fGeometryCodeSurf_micro;
    int fGeometryCode_displ_int,fGeometryCode_micro_int, fGeometryCodeSurf_displ_int;

    /** number of integration points */
    int        fNumIP_displ, fNumIPSurf_displ, fNumIP_micro, fNumIPSurf_micro;
    int knum_d_state, knum_i_state, knumstress, knumstrain;//,knumdispl;
    int num_sidesets;

    /*@}*/

    /** \name element displacements in local ordering */
    /*@{*/
    LocalArrayT u;                //solid displacement
    LocalArrayT u_n;         //solid displacement from time t_n
    LocalArrayT u_dot;         //solid velocity
    LocalArrayT u_dot_n;         //solid velocity from time t_n
    LocalArrayT u_dotdot;         //solid acceleration
    LocalArrayT u_dotdot_n;         //solid acceleration from time t_n
    LocalArrayT del_u;        //displacement increment including Newton-R update, i.e. del_u = u_{n+1}^{k+1} - u_n
    LocalArrayT Phi;        //micro-displacement-gradient dof
    LocalArrayT Phi_dot;        //micro-displacement-gradient first time derivative
    LocalArrayT Phi_dot_n;        //micro-displacement-gradient first time derivative from time t_n
    LocalArrayT Phi_dotdot;        //micro-displacement-gradient second derivative
    LocalArrayT Phi_dotdot_n;        //micro-displacement-gradient second derivative from time t_n
    LocalArrayT Phi_n;        //micro-displacement-gradient from time t_n
    LocalArrayT del_Phi;        //micro-displacement-gradient increment
    dArrayT                del_u_vec;          // vector form
    dArrayT                del_Phi_vec;        // vector form
    dArrayT                u_vec;          // solid displacement in vector form
    dArrayT                u_dot_vec;          // solid velocity in vector form
    dArrayT                u_dotdot_vec;          // solid acceleration in vector form
    dArrayT                Phi_vec;        // micro-displacement-gradient in vector form
    dArrayT                Phi_dot_vec;        // first derivative of micro-displacement-gradient in vector form
    dArrayT                Phi_dotdot_vec;        // second derivative of micro-displacement-gradient in vector form

    /*@}*/

    // problem size definitions
    int n_en_displ, n_en_displ_x_n_sd, n_sd_x_n_sd,n_sd_x_n_sd_x_n_sd,n_en_micro_x_n_sd_x_n_sd, n_en_micro_x_n_sd;
    int n_el, n_sd, n_sd_surf, n_en_surf;
    int n_en_micro, ndof_per_nd_micro, n_en_micro_x_ndof_per_nd_micro, ndof_per_nd_micro_x_n_sd;
    int step_number;
    int iConstitutiveModelType;
    //double Alpha;

    //name of output vector
    StringT output;

    dArrayT fForces_at_Node;
    bool bStep_Complete;

    double time, fRho, fRho_0;

    void Get_Fd_ext ( dArrayT &fFd_ext );

    //-- Material Parameters
    dArrayT fMaterial_Params;

    //-- Newmark Time Integration Parameters
    dArrayT fIntegration_Params;

    /** \name shape functions wrt to current coordinates */
    /*@{*/
    /** shape functions and derivatives. The derivatives are wrt to the
      * reference coordinates */
    ShapeFunctionT* fShapes_displ;
    ShapeFunctionT* fShapes_micro;

    /** reference coordinates */
    LocalArrayT fInitCoords_displ, fInitCoords_micro;
    /** current coordinates */
    LocalArrayT fCurrCoords_displ, fCurrCoords_micro;
    /*@}*/

    /* Data Storage */
    ElementMatrixT fKdd, fKdphi;
    ElementMatrixT fKphid, fKphiphi;
    dArrayT         fFd_int;
    dArrayT         fFd_ext;
    dArrayT         fFphi_int;
    dArrayT         fFphi_ext;

    dArrayT         fGrad_disp_vector;
    dMatrixT        fGrad_disp_matrix;

    dArrayT         fTemp_nine_values;
    dArrayT         fTemp_six_values;

    dMatrixT        fDeformation_Gradient;
    dMatrixT        fRight_Cauchy_Green_tensor;
    dMatrixT        fRight_Cauchy_Green_tensor_Inverse;
    dMatrixT        fRight_Elastic_Cauchy_Green_tensor;
    dMatrixT        fRight_Elastic_Cauchy_Green_tensor_tr;

    dMatrixT        fLeft_Cauchy_Green_tensor;
    dMatrixT        fLeft_Cauchy_Green_tensor_Inverse;
    dMatrixT        fDeformation_Gradient_Inverse;
    dMatrixT        fDeformation_Gradient_Transpose;
    dMatrixT        fDefGradInv_Grad_grad;

    dMatrixT        fIdentity_matrix;
    dMatrixT        fTemp_matrix_nsd_x_nsd;
    dMatrixT        fTemp_matrix_nsd_x_nsd2;
    dMatrixT        fShapeMicro_row_matrix;


    dMatrixT        u_dotdot_column_matrix;
    dMatrixT        u_dot_column_matrix;
    dMatrixT        u_dot_column_matrix_Transpose;



    //////////////////////////////////////////////////////////
    /////DEFINITIONS FOR MICROMORPHIC MATRICES////////////////
    //////////////////////////////////////////////////////////
    double KrDelta[3][3];
    double trdeltad;
    double trdeltaEp;

    double Counter;
    dArray2DT Counter_IPs_el_n;
    dArray2DT Counter_IPs_el;
    dArrayT Counter_IPs;
    dMatrixT fIota_temp_matrix;

    //Varitional Matrices coming from the Balance of linear Momentum
    dMatrixT Tsigma_1;
    dMatrixT fG1_1;
    double SigN[3][3];//unsymmetric Cauchy stress tensor found at previous step
    dMatrixT SigN_m;
    dMatrixT Fn_m;
    dMatrixT Finv_m;
    dMatrixT deltaL;
    dMatrixT deltaL_Tr;
    dMatrixT tempSig;
    dMatrixT deltad;
    dMatrixT SigN_ar;


    dMatrixT Sigma8;
    dArray2DT SigN_IPs;
    dArray2DT Sig_IPs;
    dArray2DT SigN_IPs_n;
    dArray2DT SigN_IPs_el;
    dArray2DT Sig_IPs_el;
    dArray2DT SigN_IPs_el_n;
    dArrayT Temp_Identity_array;


    dMatrixT SPiolaN;
    dMatrixT SPiola;

    dArray2DT SPiolaN_IPs;
    dArray2DT SPiola_IPs;
    dArray2DT SPiolaN_IPs_el;
    dArray2DT SPiola_IPs_el;
    dArray2DT SPiolaN_IPs_el_n;



    dMatrixT Sigma; // unsymetric Cauchy stress tensor at current step
    double Fn[3][3];
    dMatrixT Fn_ar;
    dArray2DT Fn_ar_IPs;
    dArray2DT F_ar_IPs;
    dArray2DT Fn_ar_IPs_el;
    dArray2DT F_ar_IPs_el;
    dArray2DT Fn_ar_IPs_el_n;

    double FnInv[3][3];
    dMatrixT FnInv_ar;
    dArray2DT FnInv_ar_IPs;
    dArray2DT FInv_ar_IPs;
//  dArray2DT FInv_ar_IPs_n;
    dArray2DT FnInv_ar_IPs_el;
    dArray2DT FInv_ar_IPs_el;
    dArray2DT FnInv_ar_IPs_el_n;

    double Finv[3][3];

    double Chi[3][3];
    dMatrixT Chi_m;
    double ChiInv[3][3];
    dMatrixT ChiInv_m;
    dMatrixT ChiN_m;
    dMatrixT deltaEp;
    dMatrixT deltaNu;
    double ChiN[3][3];
    dMatrixT ChinN_m;
    dMatrixT ChiN_ar;
    dMatrixT Chi_ar;
    dArray2DT ChiN_ar_IPs;
    dArray2DT Chi_ar_IPs;
    dArray2DT ChiN_ar_IPs_n;
    dArray2DT Chi_ar_IPs_el;
    dArray2DT ChiN_ar_IPs_el;
    dArray2DT ChiN_ar_IPs_el_n;

//  double  Trial[6];
    dMatrixT fIota_w_temp_matrix;
    dMatrixT fTemp_matrix_nudof_x_nchidof;
    dMatrixT fTemp_matrix_nchidof_x_nchidof;
    dMatrixT fTemp_matrix_nchidof_x_nudof;
    dMatrixT fTemp_matrix_nudof_x_nudof;

    dMatrixT fH1_Etagrad;
    dMatrixT TransShapeDisplGrad;
    dMatrixT Var_F;
    dMatrixT GRAD_Nuw;
    dMatrixT GRAD_Nuw_Tr;
    dMatrixT Finv_w; // to create Iota_w which is different than Iota because sequence of the components in wk,l
    dMatrixT Tsigma_2;
    dMatrixT fG1_2;//not being calculated yet
    dMatrixT Tsigma_3;
    dMatrixT fG1_3;// not being calculated yet
    dMatrixT TFn_1;// to be multiplied by (lamda+Tau)
    dMatrixT fG1_4;
    dMatrixT TFn_2;// to be multiplied by (Mu+sigma)
    dMatrixT fG1_5;
    dMatrixT TFn_3;// to be multiplied by (Mu+sigma)
    dMatrixT fG1_6;
    dMatrixT TChi_1;// to be multiplied by eta
    dMatrixT fG1_7;
    dMatrixT TFn_4;// to be multiplied by eta
    dMatrixT fG1_8;
    dMatrixT TChi_2;// to be multiplied by kappa
    dMatrixT fG1_9;
    dMatrixT TFn_5;// to be multiplied by kappa
    dMatrixT fG1_10;
    dMatrixT TChi_3;// to be multiplied by nu
    dMatrixT fG1_11;
    dMatrixT TFn_6;// to be multiplied by nu
    dMatrixT fG1_12;
    dMatrixT SigCurr;
    dMatrixT fG1_13;
    dMatrixT fG1_14;
 // Variational Matrices coming from the Balance of First Moment of Momentum
    double mn[3][3][3];
    dArrayT mn_ar;
    dArray2DT mn_IPs;
    dArray2DT mn_IPs_n;
    dArray2DT mn_IPs_el;
    dArray2DT mn_IPs_el_n;

    double Mnplus1[3][3][3];
    double Gamma[3][3][3];
    double GammaN[3][3][3];
    dArrayT GammaN_ar;

    dArray2DT GammaN_IPs;
    dArray2DT GammaN_IPs_n;
    dArray2DT GammaN_IPs_el;
    dArray2DT GammaN_IPs_el_n;

    double CCof[3][3][3][3][3][3];
    double GRAD_ChiN[3][3][3];// GRADIENT  in reference configuration!
    dArrayT GRAD_Chi_ar;
    dArrayT GRAD_ChiN_ar;
    dArray2DT GRAD_ChiN_ar_IPs;
    dArray2DT GRAD_Chi_ar_IPs;
    dArray2DT GRAD_ChiN_ar_IPs_n;
    dArray2DT GRAD_ChiN_ar_IPs_el;
    dArray2DT GRAD_Chi_ar_IPs_el;
    dArray2DT GRAD_ChiN_ar_IPs_el_n;

    double GRAD_Chi[3][3][3];// GRADIENT in reference configuration!



    dMatrixT fIota_eta_temp_matrix;
    dMatrixT Finv_eta; // to create Iota_eta


    dMatrixT Etagrad;
    dMatrixT Mm_1;
    dMatrixT Mm_2;
    dMatrixT Mm_3;
    dMatrixT Mm_4;
    dMatrixT Mm_5;
    dMatrixT Mm_6;
    dMatrixT Mm_7;
    dMatrixT Mm_71;
    dMatrixT Mm_72;
    dMatrixT Mm_73;
    dMatrixT Mm_74;
    dMatrixT Mm_75;
    dMatrixT Mm_76;
    dMatrixT Mm_77;
    dMatrixT Mm_78;
    dMatrixT Mm_8;
    dMatrixT Mm_9;
    dMatrixT Mm_10;
    dMatrixT Mm_11;
    dMatrixT Mm_12;
    dMatrixT Mm_13;
    dMatrixT Mm_14;
    dMatrixT Ru_1;//u
    dMatrixT Ru_2;
    dMatrixT Ru_3;
    dMatrixT RChi_1;
    dMatrixT Ru_4;
    dMatrixT RChi_2;
    dMatrixT Ru_5;
    dMatrixT Ru_6;
    dMatrixT Ru_7;
    dMatrixT Ru_8;
    dMatrixT Ru_9;
    dMatrixT RChi_3;
    dMatrixT Rs_sigma;
    dMatrixT R_gradu;
    dMatrixT R_Capital_Gamma_Chi;
    dMatrixT CapitalLambda;

    dMatrixT sn_sigman;
    dArray2DT sn_sigman_IPs;
    dArray2DT sn_sigman_IPs_n;
    dArray2DT sn_sigman_IPs_el;
    dArray2DT sn_sigman_IPs_el_n;

    dMatrixT s_sigma;

    dMatrixT fShapeDispl_Tr;


    dMatrixT NCHI;
    dMatrixT NCHI_Tr;
//    dMatrixT NCHI_eta;//  same with the one above no need!
//    dMatrixT GRAD_NCHI_Phi;//no need for this same with GRAD_NCHI
    dMatrixT GRAD_NCHI;


    dMatrixT fH1_1;
    dMatrixT fH1_2;
    dMatrixT fH1_3;
    dMatrixT fH1_4;
    dMatrixT fH1_5;
    dMatrixT fH1_6;
    dMatrixT fH1_7;
    dMatrixT fH1_71;
    dMatrixT fH1_72;
    dMatrixT fH1_73;
    dMatrixT fH1_74;
    dMatrixT fH1_75;
    dMatrixT fH1_76;
    dMatrixT fH1_77;
    dMatrixT fH1_78;
    dMatrixT fH1_8;
    dMatrixT fH1_9;
    dMatrixT fH1_10;
    dMatrixT fH1_11;
    dMatrixT fH1_12;
    dMatrixT fH1_13;
    dMatrixT fH1_14;

    dMatrixT fH2_1;
    dMatrixT fH2_2;
    dMatrixT fH2_3;
    dMatrixT fH2_4;
    dMatrixT fH2_5;
    dMatrixT fH2_6;
    dMatrixT fH2_7;
    dMatrixT fH2_8;
    dMatrixT fH2_9;
    dMatrixT fH2_10;
    dMatrixT fH2_11;
    dMatrixT fH2_12;
    dMatrixT fH2_13;
   // dMatrixT fH2_14;

    dMatrixT fH3_1;

    dArrayT Chi_vec;
    dArrayT GRAD_Chi_vec;
    //internal Force Matrices
    dArrayT G1;
    dArrayT Uint_1;
    dArrayT Uint_1_temp;
    dArrayT Pint_1_temp;
    dArrayT Pint_2_temp;
    dArrayT Pint_3_temp;
    dArrayT Uext_1;
    dArrayT Text;
    dArrayT Gext;

    dArrayT H1;
    dArrayT Pint_1;
    dArrayT H2;
    dArrayT Pint_2;
    dArrayT H3;
    dArrayT Pint_3;
    dArrayT Hext;
    dArrayT Pext;

    dMatrixT Lambda;
    dMatrixT Omega;
    //////////////////////////////////////////////////////////
    //////FINITE STRAIN ELASTICITY MATRICES START HERE/////////
    //////////////////////////////////////////////////////////
    //dMatrixT V_1;
    dMatrixT Jmat;
    dMatrixT KJmat;
    dMatrixT SPK;
    dMatrixT devSPK;
    dMatrixT KirchhoffST;// The second Piola-Kirchhoff Matrix
    dMatrixT Temp_SPK;//temporary Matrix used in calculation of SPK
  //  dMatrixT FSF;
  //  dMatrixT LST;//Lagrangian strain tensor used in some functions to get rid of long name
    dMatrixT LagrangianStn;
    dMatrixT Elastic_LagrangianStn;
    dMatrixT Elastic_LagrangianStn_tr;

    dMatrixT MicroStnTensor;//Micro-strain tensor
    dMatrixT PSI;//deformation measure PSI=Transpose(F).chi
    dMatrixT ChiM; //Micro-deformation tensor Chi ( used a different tensor this time )
    dMatrixT I1_1;
    dMatrixT I1_2;
    dMatrixT I1_3;
    dMatrixT I1_4;
    dMatrixT I1_5;
    dMatrixT I1_6;
    dMatrixT I1_7;
    dMatrixT I2_1;
    dMatrixT I1_8;
    dMatrixT I2_2;
    dMatrixT I1_9;
    dMatrixT I2_3;
/*********************/
    dMatrixT fFJ;
    dMatrixT fJF;
    dMatrixT fJ1_1;
    dMatrixT fJ1_2;
    dMatrixT fJ1_3;
    dMatrixT fJ1_4;
    dMatrixT fJ2_1;
    dMatrixT fJ1_5;
    dMatrixT fJ2_2;
    dMatrixT fJ1_6;
    dMatrixT fJ2_3;

    dArrayT Vint_1;
    dArrayT Vint_1_temp;
    dArrayT Vint_2;
    dArrayT Vint_2_temp;
    dArrayT Vint_3_temp;
    dArrayT Vint_3;
    dArrayT fV1;
    dArrayT fV2;
    dArrayT fV3;
    dMatrixT fKu_1;
    dMatrixT fKu_2;
    dMatrixT fKu_3;
    dMatrixT fKu_4;
    dMatrixT fKu_5;
    dMatrixT fKu_6;
    dMatrixT fKu_7;
    dMatrixT fKuphi_1;
    dMatrixT fKu_8;
    dMatrixT fKuphi_2;
    dMatrixT fKu_9;
    dMatrixT fKuphi_3;
    dMatrixT SIGMA_S;
    dMatrixT devSIGMA_S;
    dMatrixT fKFJu;
    dMatrixT fKJFu;
    dMatrixT fKphiu_1;
    dMatrixT fKphiu_2;
    dMatrixT fKphiu_3;
    dMatrixT fKphiu_4;
    dMatrixT fKphiphi_1;
    dMatrixT fKphiu_5;
    dMatrixT fKphiphi_2;
    dMatrixT fKphiu_6;
    dMatrixT fKphiphi_3;

    dMatrixT fFM;
    dMatrixT fMF;
    dMatrixT fMchi;
    dMatrixT fEtaM;
    dMatrixT fMpu_1;
    dMatrixT fMpp_1;
    dMatrixT fMpu_2;
    dMatrixT fMpp_2;
    dMatrixT fMpu_3;
    dMatrixT fMpp_3;
    dMatrixT fMpu_4;
    dMatrixT fMpp_4;
    dMatrixT fMpu_6;
    dMatrixT fMpp_6;
    dMatrixT fKMphiu_3;
    dMatrixT fKMphiphi_3;
    dMatrixT fKMphiu_4;
    dMatrixT fKMphiphi_4;
    dMatrixT fKMphiu_6;
    dMatrixT fKMphiphi_6;

    dTensor3DT fMKLM;
    dTensor3DT fdevMKLM;
    dTensor3DT GAMMA;
    dTensor3DT GRAD_CHIM;
    dTensor3DT fTemp_tensor_n_sd_x_n_sd_x_nsd;



    dMatrixT fKEtaM;
    dMatrixT fKMFphiu;
    dMatrixT fKMchiphiphi;
    dMatrixT fKMphiu_1;
    dMatrixT fKMphiphi_1;
    dMatrixT fKMphiu_2;
    dMatrixT fKMphiphi_2;
    /* Plasticity Matrices*/
    dMatrixT PHIMATRIX;
    dTensor3DT GPHIMATRIX;

    dArrayT Vintp_1;
    dArrayT Vintp_1_temp;
    dArrayT fV1p;
    dArrayT Vintp_2;
    dArrayT Vintp_2_temp;
    dArrayT fV2p;



    dMatrixT fFp;
    dMatrixT fFp_inverse;
    dMatrixT fFp_n;
    dMatrixT fFp_n_inverse;

    dMatrixT fChip;
    dMatrixT fChip_inverse;
    dMatrixT fChip_n;
    dMatrixT fChip_n_inverse;

    dMatrixT fFe;
    dMatrixT fFe_tr;
    dMatrixT fChie;
    dMatrixT fChie_tr;
    dMatrixT fdGdS_n;
    dMatrixT fdGdS;
    dMatrixT fdFYdS;
    dMatrixT fdFYdS_n;
    
    /***************************************************/
    /*****Micro-scale plasticity matrices **************/  
    /***************************************************/ 
    dMatrixT fdFYchidSIGMA_S;       
    dMatrixT fdGchidSIGMA_S;
    dMatrixT fdGchidSIGMA_S_n_transpose;
    dMatrixT fdGchidSIGMA_S_n;
    dMatrixT PSIe_n_inverse;    
    dMatrixT PSIe_n;    
    dMatrixT PSIe_inverse;   
    dMatrixT PSIe_tr;
    dMatrixT PSIe;
    dMatrixT fCchie;    
    dMatrixT fCchie_n;
    dMatrixT fCchie_n_inverse;    
    dMatrixT dChipdDgammachi;
    dMatrixT dChiedDgammachi;
    dMatrixT dEpsilonedDelgammachi;
    dMatrixT dSdDelgammachi;
    dMatrixT ddevSdDelgammachi;    
    dMatrixT dSIGMA_SdDelgamma;
    dMatrixT dSIGMA_SdDelgammachi;
    dMatrixT ddevSIGMA_SdDelgamma;    
    dMatrixT ddevSIGMA_SdDelgammachi;
    
    
    dMatrixT LocalConsistentTangent;
    dMatrixT LocalConsistentTangentInverse;    
    dArrayT LocalRHSVector;
    dArrayT fdelDelgammaVector;    
    
    
    /***************************************************/
    /***************************************************/
    /***************************************************/         
    
    

    dMatrixT fRight_Cauchy_Green_tensor_tr;
    dMatrixT fLagrangian_strain_tensor_tr;
    dMatrixT Elastic_MicroStnTensor_tr;
    dMatrixT Elastic_MicroStnTensor;
    dMatrixT fMicroRight_Cauchy_Green_tensor;        
    dMatrixT fMicroRight_Cauchy_Green_tensor_tr;    
    dMatrixT fMicroRight_Elastic_Cauchy_Green_tensor;        
    dMatrixT fMicroRight_Elastic_Cauchy_Green_tensor_tr;    
    //dMatrixT fMicroRight_Cauchy_Green_tensor;    
    //dMatrixT fMicroRight_Cauchy_Green_tensor_tr;

    dMatrixT fSPK_tr;
    dMatrixT fdevSPK_tr;
    //dMatrixT fSPK;
    dMatrixT dSdDelgamma;
    dMatrixT ddevSdDelgamma;
    dMatrixT dEedDelgamma;
    dMatrixT dEpsilonedDelgamma;

    //dMatrixT fSIGMA_S;
    dMatrixT SIGMA_S_tr;
    dMatrixT devSIGMA_S_tr;

    dMatrixT fFeT;


    dMatrixT fA1;
    dMatrixT fA2;

    dMatrixT fN1;
    dMatrixT fD1;

    double fYield_function,fYield_function_tr,dFYdc;
    double fMicroYield_function,fMicroYield_function_tr;
    double devfSPKinv,devfSPKinv_tr;
    double devSIGMA_S_inv,devSIGMA_S_inv_tr;    
    double fDelgamma, fdelDelgamma,dFYdDelgamma,dFYdDelgammachi;
    double fDelgammachi, fdelDelgammachi,dFYchidDelgammachi,dFYchidDelgamma;
    double dPdDelgamma,dcdDelgamma,Temp_inv,press,dinvSdDelgamma;
    double dPdDelgammachi,ddevSdDelgamma_inv,ddevSdDelgammachi_inv;
    double dPchidDelgammachi,dcchidDelgammachi,ddevSIGMA_SdDelgammachi_inv;
    double dPchidDelgamma,ddevSIGMA_SdDelgamma_inv;
    int iter_count, global_iteration;
    double Aphi,Bphi,Apsi,Bpsi;
    double Aphi_chi,Bphi_chi,Apsi_chi,Bpsi_chi;
    double Beta;

    /* for local Newton-Raphson iteration */
   int iIterationMax;
   double dRelTol, dAbsTol,fConst1;

   int PlasticityCondition;

   /* for local trial yield check */
   double dYieldTrialTol;
     /* some scalars used in calculations */
   double dFYdScol1,fdFYdS_fA1,fdFYdS_fA1T,trfA1,fdFYdc,dFYdc_delc;
   double trfN1,fdFYdS_fN1,fdFYdS_fN1T;
   double fdFYchidSIGMA_S_fD1T,fdFYchidSIGMA_S_fD1,dFYchidSIGMA_Scol1,trfD1;
   double fdFYchidSIGMA_S_fA1,fdFYchidSIGMA_S_fA1T,fdFYchidSIGMA_S_fN1,fdFYchidSIGMA_S_fN1T;
   double fConst4,fdFYchidcchi;
   double mean_stress_tr,mean_stress;
   
   
   /* for coupled solution*/  
   double fdFYdS_fD1,fdFYdS_fD1T;
   double fConst2,fConst3;
   double Comp11,Comp12,Comp21,Comp22;
   double Comp1,Comp2,Comp3,Comp4;   

    dMatrixT IJp_1;
    dMatrixT fKu_IJp_1;
    dMatrixT IJp_2;
    dMatrixT fKu_IJp_2;
    dMatrixT IJp_3;
    dMatrixT fKu_IJp_3;
    dMatrixT IJp_4;
    dMatrixT fKu_IJp_4;
    dMatrixT IJp_5;
    dMatrixT fKu_IJp_5;
    dMatrixT IJp_6;
    dMatrixT fKu_IJp_6;
    dMatrixT IJp_7;
    dMatrixT fKuphi_IJp_7;
    dMatrixT IJp_8;
    dMatrixT fKuphi_IJp_8;
    dMatrixT IJp_9;
    dMatrixT fKuphi_IJp_9;
    /* Matrices from coupling */
    dMatrixT IJp_10;
    dMatrixT fKu_IJp_10;
    dMatrixT IJp_11;
    dMatrixT fKu_IJp_11;
    dMatrixT IJp_12;
    dMatrixT fKu_IJp_12;
    dMatrixT IJp_13;
    dMatrixT fKuphi_IJp_13;    
    dMatrixT IJp_14;
    dMatrixT fKu_IJp_14;
    dMatrixT IJp_15;
    dMatrixT fKuphi_IJp_15;       
    dMatrixT IJp_16;
    dMatrixT fKu_IJp_16; 
    dMatrixT IJp_17;
    dMatrixT fKuphi_IJp_17;                         

    dMatrixT I1e_1;
    dMatrixT fKu_I1e_1;

    dMatrixT I2e_1;
    dMatrixT fKu_I2e_1;
    dMatrixT I2p_2;
    dMatrixT fKu_I2p_2;
    dMatrixT I2p_3;
    dMatrixT fKu_I2p_3;
    dMatrixT I2p_4;
    dMatrixT fKu_I2p_4;
    dMatrixT I2p_5;
    dMatrixT fKu_I2p_5;
    dMatrixT I2p_6;
    dMatrixT fKu_I2p_6;
    dMatrixT I2p_7;
    dMatrixT fKu_I2p_7;
    dMatrixT I2p_8;
    dMatrixT fKuphi_I2p_8;
    dMatrixT I2p_9;
    dMatrixT fKuphi_I2p_9;
    dMatrixT I2p_10;
    dMatrixT fKuphi_I2p_10;
    /* Matrices from couplinb */
    dMatrixT I2p_11;
    dMatrixT fKu_I2p_11;
    dMatrixT I2p_12;
    dMatrixT fKu_I2p_12;
    dMatrixT I2p_13;
    dMatrixT fKu_I2p_13;
    dMatrixT I2p_14;
    dMatrixT fKuphi_I2p_14;    
    dMatrixT I2p_15;
    dMatrixT fKu_I2p_15;
    dMatrixT I2p_16;
    dMatrixT fKuphi_I2p_16;        
    dMatrixT I2p_17;
    dMatrixT fKu_I2p_17;    
    dMatrixT I2p_18;
    dMatrixT fKuphi_I2p_18;                            
    
    
    
    
    

    dMatrixT I3e_1;
    dMatrixT fKu_I3e_1;
    dMatrixT I3e_2;
    dMatrixT fKu_I3e_2;
    dMatrixT I3e_3;
    dMatrixT fKu_I3e_3;
    dMatrixT I3p_4;
    dMatrixT fKu_I3p_4;
    dMatrixT I3p_5;
    dMatrixT fKu_I3p_5;
    dMatrixT I3p_6;
    dMatrixT fKu_I3p_6;
    dMatrixT I3p_7;
    dMatrixT fKu_I3p_7;
    dMatrixT I3p_8;
    dMatrixT fKu_I3p_8;
    dMatrixT I3p_9;
    dMatrixT fKu_I3p_9;
    dMatrixT I3p_10;
    dMatrixT fKu_I3p_10;
    dMatrixT I3p_11;
    dMatrixT fKu_I3p_11;
    dMatrixT I3p_12;
    dMatrixT fKu_I3p_12;
    dMatrixT I3e_13;
    dMatrixT fKu_I3e_13;
    dMatrixT I3e_14;
    dMatrixT fKu_I3e_14;
    dMatrixT I3e_15;
    dMatrixT fKu_I3e_15;
    dMatrixT I3p_16;
    dMatrixT fKu_I3p_16;
    dMatrixT I3p_17;
    dMatrixT fKu_I3p_17;
    dMatrixT I3p_18;
    dMatrixT fKu_I3p_18;
    dMatrixT I3p_19;
    dMatrixT fKu_I3p_19;
    dMatrixT I3p_20;
    dMatrixT fKu_I3p_20;
    dMatrixT I3p_21;
    dMatrixT fKu_I3p_21;
    dMatrixT I3p_22;
    dMatrixT fKu_I3p_22;
    dMatrixT I3p_23;
    dMatrixT fKu_I3p_23;
    dMatrixT I3p_24;
    dMatrixT fKu_I3p_24;
    dMatrixT I3p_25;
    dMatrixT fKu_I3p_25;
    dMatrixT I3p_26;
    dMatrixT fKu_I3p_26;
    dMatrixT I3p_27;
    dMatrixT fKu_I3p_27;
    dMatrixT I3p_28;
    dMatrixT fKu_I3p_28;
    dMatrixT I3p_29;
    dMatrixT fKu_I3p_29;
    dMatrixT I3p_30;
    dMatrixT fKu_I3p_30;
    dMatrixT I3p_31;
    dMatrixT fKu_I3p_31;
    dMatrixT I3p_32;
    dMatrixT fKu_I3p_32;
    dMatrixT I3p_33;
    dMatrixT fKu_I3p_33;
    dMatrixT I3p_34;
    dMatrixT fKu_I3p_34;
    dMatrixT I3p_35;
    dMatrixT fKu_I3p_35;
    dMatrixT I3p_36;
    dMatrixT fKu_I3p_36;
    dMatrixT I3e_37;
    dMatrixT fKuphi_I3e_37;
    dMatrixT I3e_38;
    dMatrixT fKuphi_I3e_38;
    dMatrixT I3e_39;
    dMatrixT fKuphi_I3e_39;
    dMatrixT I3p_40;
    dMatrixT fKuphi_I3p_40;
    dMatrixT I3p_41;
    dMatrixT fKuphi_I3p_41;
    dMatrixT I3p_42;
    dMatrixT fKuphi_I3p_42;
    dMatrixT I3p_43;
    dMatrixT fKuphi_I3p_43;
    dMatrixT I3p_44;
    dMatrixT fKuphi_I3p_44;
    dMatrixT I3p_45;
    dMatrixT fKuphi_I3p_45;
    dMatrixT I3p_46;
    dMatrixT fKuphi_I3p_46;
    dMatrixT I3p_47;
    dMatrixT fKuphi_I3p_47;
    dMatrixT I3p_48;
    dMatrixT fKuphi_I3p_48;
    dMatrixT I3p_49;
    dMatrixT fKuphi_I3p_49;
    dMatrixT I3p_50;
    dMatrixT fKuphi_I3p_50;
    dMatrixT I3p_51;
    dMatrixT fKuphi_I3p_51;
    dMatrixT I3p_52;
    dMatrixT fKuphi_I3p_52;
    dMatrixT I3p_53;
    dMatrixT fKuphi_I3p_53;
    dMatrixT I3p_54;
    dMatrixT fKuphi_I3p_54;
    /* Matrices from Del(delgammachi) */
    dMatrixT I3p_55;
    dMatrixT fKu_I3p_55;  
    dMatrixT I3p_56;
    dMatrixT fKu_I3p_56;     
    dMatrixT I3p_57;
    dMatrixT fKu_I3p_57; 
    dMatrixT I3p_58;
    dMatrixT fKuphi_I3p_58;
    dMatrixT I3p_59;
    dMatrixT fKu_I3p_59;     
    dMatrixT I3p_60;
    dMatrixT fKuphi_I3p_60;
    dMatrixT I3p_61;
    dMatrixT fKu_I3p_61;     
    dMatrixT I3p_62;
    dMatrixT fKuphi_I3p_62;            
    dMatrixT I3p_63;
    dMatrixT fKu_I3p_63;  
    dMatrixT I3p_64;
    dMatrixT fKu_I3p_64;     
    dMatrixT I3p_65;
    dMatrixT fKu_I3p_65; 
    dMatrixT I3p_66;
    dMatrixT fKuphi_I3p_66;
    dMatrixT I3p_67;
    dMatrixT fKu_I3p_67;     
    dMatrixT I3p_68;
    dMatrixT fKuphi_I3p_68;
    dMatrixT I3p_69;
    dMatrixT fKu_I3p_69;     
    dMatrixT I3p_70;
    dMatrixT fKuphi_I3p_70;      
    dMatrixT I3p_71;
    dMatrixT fKu_I3p_71;  
    dMatrixT I3p_72;
    dMatrixT fKu_I3p_72;     
    dMatrixT I3p_73;
    dMatrixT fKu_I3p_73; 
    dMatrixT I3p_74;
    dMatrixT fKuphi_I3p_74;
    dMatrixT I3p_75;
    dMatrixT fKu_I3p_75;     
    dMatrixT I3p_76;
    dMatrixT fKuphi_I3p_76;
    dMatrixT I3p_77;
    dMatrixT fKu_I3p_77;     
    dMatrixT I3p_78;
    dMatrixT fKuphi_I3p_78;    
    /* Matrices from coupling*/
    dMatrixT I3p_79;   
    dMatrixT fKu_I3p_79;   
    dMatrixT I3p_80; 
    dMatrixT fKu_I3p_80;             
    dMatrixT I3p_81;   
    dMatrixT fKu_I3p_81;       
    dMatrixT I3p_82;  
    dMatrixT fKuphi_I3p_82;                
    dMatrixT I3p_83;
    dMatrixT fKu_I3p_83;        
    dMatrixT I3p_84; 
    dMatrixT fKuphi_I3p_84;                                 
    dMatrixT I3p_85; 
    dMatrixT fKu_I3p_85;                               
    dMatrixT I3p_86;     
    dMatrixT fKuphi_I3p_86;      
    dMatrixT I3p_87;   
    dMatrixT fKu_I3p_87;   
    dMatrixT I3p_88; 
    dMatrixT fKu_I3p_88;             
    dMatrixT I3p_89;   
    dMatrixT fKu_I3p_89;       
    dMatrixT I3p_90;  
    dMatrixT fKuphi_I3p_90;                
    dMatrixT I3p_91;
    dMatrixT fKu_I3p_91;        
    dMatrixT I3p_92; 
    dMatrixT fKuphi_I3p_92;                                 
    dMatrixT I3p_93; 
    dMatrixT fKu_I3p_93;                               
    dMatrixT I3p_94;     
    dMatrixT fKuphi_I3p_94; 
    dMatrixT I3p_95;   
    dMatrixT fKu_I3p_95;   
    dMatrixT I3p_96; 
    dMatrixT fKu_I3p_96;             
    dMatrixT I3p_97;   
    dMatrixT fKu_I3p_97;       
    dMatrixT I3p_98;  
    dMatrixT fKuphi_I3p_98;                
    dMatrixT I3p_99;
    dMatrixT fKu_I3p_99;        
    dMatrixT I3p_100; 
    dMatrixT fKuphi_I3p_100;                                 
    dMatrixT I3p_101; 
    dMatrixT fKu_I3p_101;                               
    dMatrixT I3p_102;     
    dMatrixT fKuphi_I3p_102;
    dMatrixT I3p_103; 
    dMatrixT fKu_I3p_103;      
    dMatrixT I3p_104; 
    dMatrixT fKu_I3p_104;
    dMatrixT I3p_105; 
    dMatrixT fKu_I3p_105;
    dMatrixT I3p_106;     
    dMatrixT fKuphi_I3p_106;  
    dMatrixT I3p_107; 
    dMatrixT fKu_I3p_107;     
    dMatrixT I3p_108;     
    dMatrixT fKuphi_I3p_108; 
    dMatrixT I3p_109; 
    dMatrixT fKu_I3p_109;     
    dMatrixT I3p_110;     
    dMatrixT fKuphi_I3p_110;           
    dMatrixT I3p_111;   
    dMatrixT fKu_I3p_111;   
    dMatrixT I3p_112; 
    dMatrixT fKu_I3p_112;             
    dMatrixT I3p_113;   
    dMatrixT fKu_I3p_113;       
    dMatrixT I3p_114;  
    dMatrixT fKuphi_I3p_114;                
    dMatrixT I3p_115;
    dMatrixT fKu_I3p_115;        
    dMatrixT I3p_116; 
    dMatrixT fKuphi_I3p_116;                                 
    dMatrixT I3p_117; 
    dMatrixT fKu_I3p_117;                               
    dMatrixT I3p_118;     
    dMatrixT fKuphi_I3p_118;  
    dMatrixT I3p_119;   
    dMatrixT fKu_I3p_119;   
    dMatrixT I3p_120; 
    dMatrixT fKu_I3p_120;             
    dMatrixT I3p_121;   
    dMatrixT fKu_I3p_121;       
    dMatrixT I3p_122;  
    dMatrixT fKuphi_I3p_122;                
    dMatrixT I3p_123;
    dMatrixT fKu_I3p_123;        
    dMatrixT I3p_124; 
    dMatrixT fKuphi_I3p_124;                                 
    dMatrixT I3p_125; 
    dMatrixT fKu_I3p_125;                               
    dMatrixT I3p_126;     
    dMatrixT fKuphi_I3p_126;  
    dMatrixT I3p_127; 
    dMatrixT fKu_I3p_127;      
    dMatrixT I3p_128; 
    dMatrixT fKu_I3p_128;
    dMatrixT I3p_129; 
    dMatrixT fKu_I3p_129;
    dMatrixT I3p_130;     
    dMatrixT fKuphi_I3p_130;  
    dMatrixT I3p_131; 
    dMatrixT fKu_I3p_131;     
    dMatrixT I3p_132;     
    dMatrixT fKuphi_I3p_132; 
    dMatrixT I3p_133; 
    dMatrixT fKu_I3p_133;     
    dMatrixT I3p_134;     
    dMatrixT fKuphi_I3p_134;     
    dMatrixT I3p_135; 
    dMatrixT fKu_I3p_135;      
    dMatrixT I3p_136; 
    dMatrixT fKu_I3p_136;
    dMatrixT I3p_137; 
    dMatrixT fKu_I3p_137;
    dMatrixT I3p_138;     
    dMatrixT fKuphi_I3p_138;  
    dMatrixT I3p_139; 
    dMatrixT fKu_I3p_139;     
    dMatrixT I3p_140;     
    dMatrixT fKuphi_I3p_140; 
    dMatrixT I3p_141; 
    dMatrixT fKu_I3p_141;     
    dMatrixT I3p_142;     
    dMatrixT fKuphi_I3p_142;                                                                       

               

    dMatrixT I4e_1;
    dMatrixT fKu_I4e_1;
    dMatrixT I4p_2;
    dMatrixT fKu_I4p_2;
    dMatrixT I4p_3;
    dMatrixT fKu_I4p_3;
    dMatrixT I4p_4;
    dMatrixT fKu_I4p_4;
    dMatrixT I4p_5;
    dMatrixT fKu_I4p_5;
    dMatrixT I4p_6;
    dMatrixT fKu_I4p_6;
    dMatrixT I4p_7;
    dMatrixT fKu_I4p_7;
    dMatrixT I4p_8;
    dMatrixT fKuphi_I4p_8;
    dMatrixT I4p_9;
    dMatrixT fKuphi_I4p_9;
    dMatrixT I4p_10;
    dMatrixT fKuphi_I4p_10;
    /* Matrices from coupling*/    
    dMatrixT I4p_11;
    dMatrixT fKu_I4p_11;   
    dMatrixT I4p_12;
    dMatrixT fKu_I4p_12;
    dMatrixT I4p_13;
    dMatrixT fKu_I4p_13;                 
    dMatrixT I4p_14;
    dMatrixT fKuphi_I4p_14;
    dMatrixT I4p_15;
    dMatrixT fKu_I4p_15;        
    dMatrixT I4p_16;
    dMatrixT fKuphi_I4p_16;
    dMatrixT I4p_17;
    dMatrixT fKu_I4p_17;        
    dMatrixT I4p_18;
    dMatrixT fKuphi_I4p_18;            

    dMatrixT IIJp_1;
    dMatrixT fKphiu_IIJp_1;
    dMatrixT IIJp_2;
    dMatrixT fKphiu_IIJp_2;
    dMatrixT IIJp_3;
    dMatrixT fKphiu_IIJp_3;
    dMatrixT IIJp_4;
    dMatrixT fKphiu_IIJp_4;
    dMatrixT IIJp_5;
    dMatrixT fKphiu_IIJp_5;
    dMatrixT IIJp_6;
    dMatrixT fKphiphi_IIJp_6;
    dMatrixT IIJp_7;
    dMatrixT fKphiphi_IIJp_7;
    dMatrixT IIJp_8;
    dMatrixT fKphiphi_IIJp_8;
    /* Matrices from coupling*/
    dMatrixT IIJp_9;
    dMatrixT fKphiu_IIJp_9;
    dMatrixT IIJp_10;
    dMatrixT fKphiu_IIJp_10;
    dMatrixT IIJp_11;
    dMatrixT fKphiu_IIJp_11;            
    dMatrixT IIJp_12;
    dMatrixT fKphiphi_IIJp_12;
    dMatrixT IIJp_13;
    dMatrixT fKphiu_IIJp_13;      
    dMatrixT IIJp_14;
    dMatrixT fKphiphi_IIJp_14;
    dMatrixT IIJp_15;
    dMatrixT fKphiu_IIJp_15;      
    dMatrixT IIJp_16;
    dMatrixT fKphiphi_IIJp_16;            

    dMatrixT II2e_1;
    dMatrixT fKphiu_II2e_1;
    dMatrixT II2p_2;
    dMatrixT fKphiu_II2p_2;
    dMatrixT II2p_3;
    dMatrixT fKphiu_II2p_3;
    dMatrixT II2p_4;
    dMatrixT fKphiu_II2p_4;
    dMatrixT II2p_5;
    dMatrixT fKphiu_II2p_5;
    dMatrixT II2p_6;
    dMatrixT fKphiu_II2p_6;
    dMatrixT II2p_7;
    dMatrixT fKphiphi_II2p_7;
    dMatrixT II2p_8;
    dMatrixT fKphiphi_II2p_8;
    dMatrixT II2p_9;
    dMatrixT fKphiphi_II2p_9;
    /* Matrices from coupling*/    
    dMatrixT II2p_10;
    dMatrixT fKphiu_II2p_10;
    dMatrixT II2p_11;
    dMatrixT fKphiu_II2p_11;
    dMatrixT II2p_12;
    dMatrixT fKphiu_II2p_12;
    dMatrixT II2p_13;
    dMatrixT fKphiphi_II2p_13;    
    dMatrixT II2p_14;
    dMatrixT fKphiu_II2p_14;    
    dMatrixT II2p_15;
    dMatrixT fKphiphi_II2p_15; 
    dMatrixT II2p_16;
    dMatrixT fKphiu_II2p_16;
    dMatrixT II2p_17;
    dMatrixT fKphiphi_II2p_17;                          



    dMatrixT II3e_1;
    dMatrixT fKphiu_II3e_1;
    dMatrixT II3e_2;
    dMatrixT fKphiu_II3e_2;
    dMatrixT II3e_3;
    dMatrixT fKphiu_II3e_3;
    dMatrixT II3e_4;
    dMatrixT fKphiu_II3e_4;
    dMatrixT II3e_5;
    dMatrixT fKphiu_II3e_5;
    dMatrixT II3p_6;
    dMatrixT fKphiu_II3p_6;
    dMatrixT II3p_7;
    dMatrixT fKphiu_II3p_7;
    dMatrixT II3p_8;
    dMatrixT fKphiu_II3p_8;
    dMatrixT II3p_9;
    dMatrixT fKphiu_II3p_9;
    dMatrixT II3p_10;
    dMatrixT fKphiu_II3p_10;
    dMatrixT II3p_11;
    dMatrixT fKphiu_II3p_11;
    dMatrixT II3p_12;
    dMatrixT fKphiu_II3p_12;
    dMatrixT II3p_13;
    dMatrixT fKphiu_II3p_13;
    dMatrixT II3p_14;
    dMatrixT fKphiu_II3p_14;
    dMatrixT II3p_15;
    dMatrixT fKphiu_II3p_15;
    dMatrixT II3p_16;
    dMatrixT fKphiu_II3p_16;
    dMatrixT II3p_17;
    dMatrixT fKphiu_II3p_17;
    dMatrixT II3p_18;
    dMatrixT fKphiu_II3p_18;
    dMatrixT II3p_19;
    dMatrixT fKphiu_II3p_19;
    dMatrixT II3p_20;
    dMatrixT fKphiu_II3p_20;
    dMatrixT II3p_21;
    dMatrixT fKphiu_II3p_21;
    dMatrixT II3p_22;
    dMatrixT fKphiu_II3p_22;
    dMatrixT II3p_23;
    dMatrixT fKphiu_II3p_23;
    dMatrixT II3p_24;
    dMatrixT fKphiu_II3p_24;
    dMatrixT II3p_25;
    dMatrixT fKphiu_II3p_25;
    dMatrixT II3p_26;
    dMatrixT fKphiu_II3p_26;
    dMatrixT II3p_27;
    dMatrixT fKphiu_II3p_27;
    dMatrixT II3p_28;
    dMatrixT fKphiu_II3p_28;
    dMatrixT II3p_29;
    dMatrixT fKphiu_II3p_29;
    dMatrixT II3p_30;
    dMatrixT fKphiu_II3p_30;
    dMatrixT II3e_31;
    dMatrixT fKphiphi_II3e_31;
    dMatrixT II3e_32;
    dMatrixT fKphiphi_II3e_32;
    dMatrixT II3e_33;
    dMatrixT fKphiphi_II3e_33;
    dMatrixT II3p_34;
    dMatrixT fKphiphi_II3p_34;
    dMatrixT II3p_35;
    dMatrixT fKphiphi_II3p_35;
    dMatrixT II3p_36;
    dMatrixT fKphiphi_II3p_36;
    dMatrixT II3p_37;
    dMatrixT fKphiphi_II3p_37;
    dMatrixT II3p_38;
    dMatrixT fKphiphi_II3p_38;
    dMatrixT II3p_39;
    dMatrixT fKphiphi_II3p_39;
    dMatrixT II3p_40;
    dMatrixT fKphiphi_II3p_40;
    dMatrixT II3p_41;
    dMatrixT fKphiphi_II3p_41;
    dMatrixT II3p_42;
    dMatrixT fKphiphi_II3p_42;
    dMatrixT II3p_43;
    dMatrixT fKphiphi_II3p_43;
    dMatrixT II3p_44;
    dMatrixT fKphiphi_II3p_44;
    dMatrixT II3p_45;
    dMatrixT fKphiphi_II3p_45;
    dMatrixT II3p_46;
    dMatrixT fKphiphi_II3p_46;
    dMatrixT II3p_47;
    dMatrixT fKphiphi_II3p_47;
    dMatrixT II3p_48;
    dMatrixT fKphiphi_II3p_48;
    /* Matrices from Del(delgammachi) */
    dMatrixT II3p_49;
    dMatrixT fKphiu_II3p_49;
    dMatrixT II3p_50;
    dMatrixT fKphiu_II3p_50;
    dMatrixT II3p_51;
    dMatrixT fKphiu_II3p_51; 
    dMatrixT II3p_52;
    dMatrixT fKphiphi_II3p_52;
    dMatrixT II3p_53;
    dMatrixT fKphiu_II3p_53; 
    dMatrixT II3p_54;
    dMatrixT fKphiphi_II3p_54;
    dMatrixT II3p_55;
    dMatrixT fKphiu_II3p_55;     
    dMatrixT II3p_56;
    dMatrixT fKphiphi_II3p_56;            
    dMatrixT II3p_57;
    dMatrixT fKphiu_II3p_57;
    dMatrixT II3p_58;
    dMatrixT fKphiu_II3p_58;
    dMatrixT II3p_59;
    dMatrixT fKphiu_II3p_59; 
    dMatrixT II3p_60;
    dMatrixT fKphiphi_II3p_60;
    dMatrixT II3p_61;
    dMatrixT fKphiu_II3p_61; 
    dMatrixT II3p_62;
    dMatrixT fKphiphi_II3p_62;
    dMatrixT II3p_63;
    dMatrixT fKphiu_II3p_63;     
    dMatrixT II3p_64;
    dMatrixT fKphiphi_II3p_64;  
    dMatrixT II3p_65;   
    dMatrixT fKphiu_II3p_65;
    dMatrixT II3p_66;
    dMatrixT fKphiu_II3p_66;
    dMatrixT II3p_67;
    dMatrixT fKphiu_II3p_67; 
    dMatrixT II3p_68;
    dMatrixT fKphiphi_II3p_68;
    dMatrixT II3p_69;
    dMatrixT fKphiu_II3p_69; 
    dMatrixT II3p_70;
    dMatrixT fKphiphi_II3p_70;
    dMatrixT II3p_71;
    dMatrixT fKphiu_II3p_71;     
    dMatrixT II3p_72;
    dMatrixT fKphiphi_II3p_72;  
    /* Matrices from coupling */
    dMatrixT II3p_73;
    dMatrixT fKphiu_II3p_73;
    dMatrixT II3p_74;
    dMatrixT fKphiu_II3p_74;
    dMatrixT II3p_75;
    dMatrixT fKphiu_II3p_75;   
    dMatrixT II3p_76;
    dMatrixT fKphiphi_II3p_76; 
    dMatrixT II3p_77;
    dMatrixT fKphiu_II3p_77; 
    dMatrixT II3p_78;
    dMatrixT fKphiphi_II3p_78;
    dMatrixT II3p_79;
    dMatrixT fKphiu_II3p_79;     
    dMatrixT II3p_80;
    dMatrixT fKphiphi_II3p_80;
    dMatrixT II3p_81;
    dMatrixT fKphiu_II3p_81;
    dMatrixT II3p_82;
    dMatrixT fKphiu_II3p_82;
    dMatrixT II3p_83;
    dMatrixT fKphiu_II3p_83;   
    dMatrixT II3p_84;
    dMatrixT fKphiphi_II3p_84; 
    dMatrixT II3p_85;
    dMatrixT fKphiu_II3p_85; 
    dMatrixT II3p_86;
    dMatrixT fKphiphi_II3p_86;
    dMatrixT II3p_87;
    dMatrixT fKphiu_II3p_87;     
    dMatrixT II3p_88;
    dMatrixT fKphiphi_II3p_88;             


           

    dMatrixT II4e_1;
    dMatrixT fKphiu_II4e_1;
    dMatrixT II4p_2;
    dMatrixT fKphiu_II4p_2;
    dMatrixT II4p_3;
    dMatrixT fKphiu_II4p_3;
    dMatrixT II4p_4;
    dMatrixT fKphiu_II4p_4;
    dMatrixT II4p_5;
    dMatrixT fKphiu_II4p_5;
    dMatrixT II4p_6;
    dMatrixT fKphiu_II4p_6;
    dMatrixT II4p_7;
    dMatrixT fKphiphi_II4p_7;
    dMatrixT II4p_8;
    dMatrixT fKphiphi_II4p_8;
    dMatrixT II4p_9;
    dMatrixT fKphiphi_II4p_9;


    dMatrixT dFedDelgamma;
    dMatrixT dFedDelgammachi;    

    dArray2DT   fState_variables_IPs;
    dArray2DT   fState_variables_Elements_IPs;
    dArray2DT   fState_variables_n_IPs;
    dArray2DT   fState_variables_n_Elements_IPs;

    dArray2DT   fFp_IPs;
    dArray2DT   fFp_Elements_IPs;
    dArray2DT   fFp_n_IPs;
    dArray2DT   fFp_n_Elements_IPs;

    dArray2DT   fCe_IPs;
    dArray2DT   fCe_Elements_IPs;

    dMatrixT    fCe_n;
    dMatrixT    fCe_n_inverse;
    dMatrixT    fdGdS_n_transpose;
    dArray2DT   fCe_n_IPs;
    dArray2DT   fCe_n_Elements_IPs;

    dArray2DT   fdGdS_IPs;
    dArray2DT   fdGdS_n_IPs;


    dArray2DT   fdGdS_Elements_IPs;
    dArray2DT   fdGdS_n_Elements_IPs;


    dArray2DT   fdFYdS_IPs;
    dArray2DT   fdFYdS_n_IPs;


    dArray2DT   fdFYdS_Elements_IPs;
    dArray2DT   fdFYdS_n_Elements_IPs;

    /* Micromorphic additions */
    dArray2DT   fdGchidSIGMA_S_IPs;
    dArray2DT   fdGchidSIGMA_S_Elements_IPs;        
    dArray2DT   fdGchidSIGMA_S_n_IPs;
    dArray2DT   fdGchidSIGMA_S_n_Elements_IPs;    
    
    dArray2DT   PSIe_IPs;
    dArray2DT   PSIe_Elements_IPs;   
    dArray2DT   PSIe_n_IPs;
    dArray2DT   PSIe_n_Elements_IPs;   
    
    dArray2DT   fChip_IPs;
    dArray2DT   fChip_n_IPs;

    dArray2DT   fChip_Elements_IPs;
    dArray2DT   fChip_n_Elements_IPs;

    dArray2DT   fCchie_IPs;
    dArray2DT   fCchie_Elements_IPs;    
    dArray2DT   fCchie_n_IPs;
    dArray2DT   fCchie_n_Elements_IPs;


/////stress invariants variables////////
 double Cauchy_inv;
 double trS;//trace of SPK
 double invdevS;// Norm of dev. part of SPK
 double Rel_stres_inv;
 double trSIGMA_S;
 double invdevSIGMA_S;// Norm of dev. part of SIGMA_S
 double Higher_orderT_inv;
 double invtrM;//||tr(M)|| trM is a vector
 double invdevMKLM;//||dev(M)||
 double temp_inv;
 double invPhi;
 double invGPhi;


 dMatrixT devsigma;
 dMatrixT devRelsts;
 dTensor3DT  devmklm;
 dMatrixT s_sigma_temp;
 dTensor3DT fmklm;

//////////////////////////////////////////////////
    double invJ;

    double trsigma;
    double trs_sigma;
    double trmklm;
    double invtrMKLM;
    dArrayT trvecmklm;
    dArrayT ftemp_u_element;
    dArrayT fState_variables;

    //////////////////////////////////////////////////////////
    /////DEFINITIONS FINISH HERE FOR MICROMORPHIC MATRICES////
    //////////////////////////////////////////////////////////

        /* to store fEulerian_strain_tensor_current_IP */
        dMatrixT    fEulerian_strain_tensor_current_IP;
        /* to store fEulerian_strain_IPs for each of the 27 IPs of each element */
        dArray2DT   fEulerian_strain_IPs;
        /* to store fCauchy_stress_tensor_current_IP */
        dMatrixT    fCauchy_stress_tensor_current_IP;
        /* to store fCauchy_stress_IPs for each of the 27 IPs of each element */
        dArray2DT   fCauchy_stress_IPs;
        ////////////////////////////////////
        dArray2DT  fE_values_IPs;
        dArray2DT  fVarepsilon_IPs;
        ////////////////////////////////////

        dArray2DT   fEulerian_strain_Elements_IPs;
        dArray2DT   fCauchy_stress_Elements_IPs;

        /////////////////////////////////////////
        dArray2DT  fE_values_Element_IPs;
        dArray2DT  fVarepsilon_Element_IPs;
        ////////////////////////////////////////

        dArray2DT fDisplacement_Element_IPs;
        dArrayT   fDisplacements_current_IPs;
        /* to store displacements for each of 27  IPs of each element
         ( probably it is better to do it only  for certain nodes but I do not know how to extract the coordinates of the nodes and calculate
         the shape functions a nodes)
         */
       // dArray2DT  fDisplacement_IPs;


    /** the solid displacement field */
    const FieldT* fDispl;

    /** the micro-displacement-gradient field */
    const FieldT* fMicro;

    /** \name state variable storage *
     * State variables are handled ABAQUS-style. For every iteration, the state
     * variables from the previous increment are passed to the element, which
     * updates the values in place. Each row in the array is the state variable
     * storage for all integration points for an element */
    /*@{*/
    dArray2DT fdState_new;
    dArray2DT fdState;

    iArray2DT fiState_new;
    iArray2DT fiState;
    /*@}*/

    /** \name connectivities */
    /*@{*/
    ArrayT<const iArray2DT*> fConnectivities_displ;
    ArrayT<const iArray2DT*> fConnectivities_micro;
    ArrayT<iArray2DT> fConnectivities_reduced;
    /*@}*/

    /** \name equations */
    /*@{*/
    ArrayT<iArray2DT> fEqnos_displ;
    ArrayT<iArray2DT> fEqnos_micro;
    /*@}*/

    /** \name element cards */
    /*@{*/
    AutoArrayT<ElementCardT> fElementCards_displ;
    AutoArrayT<ElementCardT> fElementCards_micro;
    /*@}*/

    /** \name output */
    /*@{*/
    /** output ID */
    int fOutputID;

    /** integration point stresses. Calculated and stored during
     * FSMicromorphic2DT::RHSDriver */
    dArray2DT fIPVariable;
    /*@}*/

    /** \name prescribed plastic gradient side set ID */
    /*@{*/
    ArrayT<StringT> fSideSetID;

    /** prescribed micro-displacement-gradient weight over the side set;
        the direction is defined by {n1,n2,n3} ?? */
    //ArrayT<iArray2DT> fMicroWght;

    /** for each side set, the global nodes on the faces in the set */
    ArrayT<iArray2DT> fMicroFaces;

    /** equation numbers for the nodes on each face */
    ArrayT<iArray2DT> fMicroFaceEqnos;

    /** side set elements */
    ArrayT<iArrayT> fSideSetElements;

    /** side set faces */
    ArrayT<iArrayT> fSideSetFaces;
    /*@}*/

    /** write output for debugging */
    /*@{*/
    /** output file stream */
    ofstreamT fs_micromorph3D_out;
    ofstreamT fs_micromorph3DMn_out;

    /** line output formating variables */
    int outputPrecision, outputFileWidth;
    /*@}*/

    void Form_solid_shape_functions(const double* &shapes_displ_X);
    void Form_Gradient_of_solid_shape_functions(const dMatrixT &fShapeDisplGrad_temp);
    void Form_micro_shape_functions(const double* &shapes_micro_X);
    void Form_deformation_gradient_tensor(void);
    void Form_Grad_grad_transformation_matrix(void);
    void Form_deformation_gradient_inv_vector(void);
    void Extract_six_values_from_symmetric_tensor(const dMatrixT &fTensor,dArrayT& fTemp_six_values);
    void Put_values_In_dArrayT_vector(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);
    void Put_values_In_Array(const dArray2DT &f2DArrayT,const int& e,const int& IP,dArrayT& fArrayT);

    //////////////////////////////////////////////////////////
    /////FUNCTIONS  FOR MICROMORPHIC MATRICES////////////////
    //////////////////////////////////////////////////////////
    //Forming the Matrices coming from the Balance of Linear Momentum
    void Form_Gamma_tensor3D(void);
    void Form_micro_deformation_tensor_Chi(void);
    void Form_Chi_inv_matrix(void);
    void Form_GRAD_Chi_matrix(void);
    void Form_double_Finv_from_Deformation_tensor_inverse(void);
    void Form_GRAD_Nuw_matrix(const dMatrixT &fShapeDisplGrad_temp);
    void Form_Finv_w_matrix(void);

    void Mapping_double_and_Array(const int condition);//
//  void Mapping_double_and_Array(double& dmat, dArrayT& fArrayT,const int& dim,const int& condition);
    void Form_deformation_tensors_arrays(const int condition);//

    void Form_KroneckerDelta_matrix(void);
    void Form_Var_F_tensor(void);
    void Form_Tsigma_1_matrix(void);
    void Form_fG1_1_matrix(void);//not defined yet
    void Form_Tsigma_2_matrix(void);
    void Form_fG1_2_matrix(void);
    void Form_Tsigma_3_matrix(void);
    void Form_fG1_3_matrix(void);
    void Form_TFn_1_matrix(void);
    void Form_fG1_4_matrix(void);
    void Form_TFn_2_matrix(void);
    void Form_fG1_5_matrix(void);
    void Form_TFn_3_matrix(void);
    void Form_fG1_6_matrix(void);
    void Form_TChi_1_matrix(void);
    void Form_fG1_7_matrix(void);
    void Form_TFn_4_matrix(void);
    void Form_fG1_8_matrix(void);
    void Form_TChi_2_matrix(void);
    void Form_fG1_9_matrix(void);
    void Form_TFn_5_matrix(void);
    void Form_fG1_10_matrix(void);
    void Form_TChi_3_matrix(void);
    void Form_fG1_11_matrix(void);
    void Form_TFn_6_matrix(void);
    void Form_fG1_12_matrix(void);
    void Form_SigCurr_matrix(void);
    void Form_fG1_13_matrix(void);
    //Forming the Matrices coming from the Balance of First Moment of Momentum
    void Form_Gradient_of_micro_shape_eta_functions(const dMatrixT &fShapeMicroGrad);
//  void Form_NCHI_eta_matrix(const dMatrixT &fShapeMicro_row_matrix);//same with the one below
    void Form_NCHI_matrix(const dMatrixT &fShapeMicro_row_matrix); //shape function matrice for micro deformatin {ETA}=[NCHI].{alpha}

    void Form_Etagrad_matrix(void);
    void Form_Mm_1_matrix(void);// needs to be multiplied by "-" and J
    void Form_Mm_2_matrix(void);// needs to be multiplied by J
    void Form_Mm_3_matrix(void);// needs to be multiplied by J
    void Form_Mm_4_matrix(void);// needs to be multiplied by J
    void Form_Mm_5_matrix(void);// needs to be multiplied by J
    void Form_Mm_6_matrix(void);// needs to be multiplied by J
    void Form_Mm_7_matrix(void);
    void Form_Mm_71_matrix(void);
    void Form_Mm_72_matrix(void);
    void Form_Mm_73_matrix(void);
    void Form_Mm_74_matrix(void);
    void Form_Mm_75_matrix(void);
    void Form_Mm_76_matrix(void);
    void Form_Mm_77_matrix(void);
    void Form_Mm_78_matrix(void);
    void Form_Mm_8_matrix(void);
    void Form_Mm_9_matrix(void);
    void Form_Mm_10_matrix(void);
    void Form_Mm_11_matrix(void);
    void Form_Mm_12_matrix(void);
    void Form_Mm_13_matrix(void);
    void Form_Mm_14_matrix(void);//something should be changed in the loop due to div(du)!!! changed done ok!
    void Form_Ru_1_matrix(void);
    void Form_Ru_2_matrix(void);
    void Form_Ru_3_matrix(void);
    void Form_RChi_1_matrix(void);// needs to be multiplied by Kappa and J
    void Form_Ru_4_matrix(void);// needs to be multiplied by Kappa and J
    void Form_RChi_2_matrix(void);//needs to be multiplied by Nu and J
    void Form_Ru_5_matrix(void);// needs to be multiplied by Nu and J
    void Form_Ru_6_matrix(void);
    void Form_Ru_7_matrix(void);
    void Form_Ru_8_matrix(void);
    void Form_Ru_9_matrix(void);
    void Form_RChi_3_matrix(void);
    void Form_Rs_sigma_matrix(void);
    void Form_R_Capital_Lambda_Chi_matrix(void);// DO NOT multiply with J !!!
    void Form_Finv_eta_matrix(void);
    void Form_CapitalLambda_matrix(void);
    void Form_CCof_tensor(void);
    void Form_SPiola_matrix(void);

//    void Form_GRAD_NCHI_Phi_matrix(const  dMatrixT &fShapeMicroGrad); no need for this

    void Form_G1_matrix(void);
    void Form_H1_matrix(void);
    void Form_H2_matrix(void);
    void Form_H3_matrix(void);


    //////////////////////////////////////////////////////////
    /////FINITE STRAIN ELASTICITY MATRICES START /////////////
    //////////////////////////////////////////////////////////
    void Form_fV1(void);
    void Form_fV2(void);
    void Form_fV3(void);
    void Form_I1_1_matrix(void);
    void Form_Second_Piola_Kirchhoff_SPK(const dMatrixT& LagStn, const dMatrixT& MicroStn);
    void Form_ChiM(void);
    void Form_I1_1(void);
    void Form_I1_2(void);
    void Form_I1_3(void);
    void Form_I1_4(void);
    void Form_I1_5(void);
    void Form_I1_6(void);
    void Form_I1_7(void);
    void Form_I2_1(void);
    void Form_I1_8(void);
    void Form_I2_2(void);
    void Form_I1_9(void);
    void Form_I2_3(void);
    void Form_SIGMA_S(void);
    void Form_fFJ(void);
    void Form_fJF(void);
    void Form_fJ1_1(void);
    void Form_fJ1_2(void);
    void Form_fJ1_3(void);
    void Form_fJ1_4(void);
    void Form_fJ2_1(void);
    void Form_fJ1_5(void);
    void Form_fJ2_2(void);
    void Form_fJ1_6(void);
    void Form_fJ2_3(void);
    void Form_fFM(void);
    void Form_fMF(void);
    void Form_fMKLM(void);
    void Form_GAMMA(void);
    void Form_fEtaM(void);
    void Form_fMchi(void);
    void Form_fMpu_1(void);
    void Form_fMpp_1(void);
    void Form_fMpu_2(void);
    void Form_fMpp_2(void);
    void Form_fMpu_3(void);
    void Form_fMpp_3(void);
    void Form_fMpu_4(void);
    void Form_fMpp_4(void);
    void Form_fMpu_6(void);
    void Form_fMpp_6(void);



    void Form_Jmat(void);


    ////////////////////////////////////////////////////////
    /////////// functions to calculate stress measures ////
    void Calculate_Cauchy_INV(void);
    void Calculate_stress_diff_INV(void);
    void Calculate_higher_order_tensor_INV(void);
    void Calculate_fmklm(void);
    void Calculate_PHI_GPHI_matrices_INV(void);
    void Caculate_invdevpart_of_Matrix(const dMatrixT &fMatrix,dMatrixT &fdevfMatrix,double devinvariant);


/* Plasticity functions */
    void  Form_fV1p(void);
    void  Form_fV2p(void);

// The first and the second terms cancel each other, we start from the thirh term which has 12 matrices
// all the matrices having "e" are the ones which direclty include the term  d(deltau)/dX
// matrices having "p" are the one which include delta(gamma)
    void Form_IJp_1(void);
    void Form_IJp_2(void);
    void Form_IJp_3(void);
    void Form_IJp_4(void);
    void Form_IJp_5(void);
    void Form_IJp_6(void);
    void Form_IJp_7(void);
    void Form_IJp_8(void);
    void Form_IJp_9(void);
    /* Matrices from coupling */
    void Form_IJp_10(void);
    void Form_IJp_11(void);
    void Form_IJp_12(void);
    void Form_IJp_13(void);
    void Form_IJp_14(void);
    void Form_IJp_15(void);
    void Form_IJp_16(void);    
    void Form_IJp_17(void);                                                        


    void Form_I1e_1(void); // the first term first matrix

    void Form_I2e_1(void); // the  second term first matrix
    void Form_I2p_2(void);
    void Form_I2p_3(void);
    void Form_I2p_4(void);
    void Form_I2p_5(void);
    void Form_I2p_6(void);
    void Form_I2p_7(void);
    void Form_I2p_8(void);
    void Form_I2p_9(void);
    void Form_I2p_10(void);
    /* Matrices from coupling */    
    void Form_I2p_11(void);
    void Form_I2p_12(void);
    void Form_I2p_13(void);
    void Form_I2p_14(void);
    void Form_I2p_15(void);
    void Form_I2p_16(void);
    void Form_I2p_17(void);
    void Form_I2p_18(void);                            



    void Form_I3e_1(void); // the third term first matrix
    void Form_I3e_2(void); //
    void Form_I3e_3(void); //
    void Form_I3p_4(void); //
    void Form_I3p_5(void); //
    void Form_I3p_6(void); //
    void Form_I3p_7(void); //
    void Form_I3p_8(void); //
    void Form_I3p_9(void); //
    void Form_I3p_10(void); //
    void Form_I3p_11(void); //
    void Form_I3p_12(void); //
    void Form_I3e_13(void); //
    void Form_I3e_14(void); //
    void Form_I3e_15(void); //
    void Form_I3p_16(void); //
    void Form_I3p_17(void); //
    void Form_I3p_18(void); //
    void Form_I3p_19(void); //
    void Form_I3p_20(void); //
    void Form_I3p_21(void); //
    void Form_I3p_22(void); //
    void Form_I3p_23(void); //
    void Form_I3p_24(void); //
    void Form_I3p_25(void); //
    void Form_I3p_26(void); //
    void Form_I3p_27(void); //
    void Form_I3p_28(void); //
    void Form_I3p_29(void); //
    void Form_I3p_30(void); //
    void Form_I3p_31(void); //
    void Form_I3p_32(void); //
    void Form_I3p_33(void); //
    void Form_I3p_34(void); //
    void Form_I3p_35(void); //
    void Form_I3p_36(void); //
    void Form_I3e_37(void); //
    void Form_I3e_38(void); //
    void Form_I3e_39(void); //
    void Form_I3p_40(void); //
    void Form_I3p_41(void); //
    void Form_I3p_42(void); //
    void Form_I3p_43(void); //
    void Form_I3p_44(void); //
    void Form_I3p_45(void); //
    void Form_I3p_46(void); //
    void Form_I3p_47(void); //
    void Form_I3p_48(void); //
    void Form_I3p_49(void); //
    void Form_I3p_50(void); //
    void Form_I3p_51(void); //
    void Form_I3p_52(void); //
    void Form_I3p_53(void); //
    void Form_I3p_54(void); //
/* Matrices from Del(delgammachi) */    
    void Form_I3p_55(void); //
    void Form_I3p_56(void); //
    void Form_I3p_57(void); //
    void Form_I3p_58(void); //
    void Form_I3p_59(void); //
    void Form_I3p_60(void); //
    void Form_I3p_61(void); //
    void Form_I3p_62(void); //
    void Form_I3p_63(void); //
    void Form_I3p_64(void); //
    void Form_I3p_65(void); //
    void Form_I3p_66(void); //
    void Form_I3p_67(void); //
    void Form_I3p_68(void); //
    void Form_I3p_69(void); //
    void Form_I3p_70(void); //     
    void Form_I3p_71(void); //
    void Form_I3p_72(void); //
    void Form_I3p_73(void); //
    void Form_I3p_74(void); //
    void Form_I3p_75(void); //
    void Form_I3p_76(void); //
    void Form_I3p_77(void); //
    void Form_I3p_78(void); //     
    /* Matrices from coupling*/
    void Form_I3p_79(void); //     
    void Form_I3p_80(void); //    
    void Form_I3p_81(void); //                                            
    void Form_I3p_82(void); //                                            
    void Form_I3p_83(void); //                                            
    void Form_I3p_84(void); //                                            
    void Form_I3p_85(void); //  
    void Form_I3p_86(void); //   
    void Form_I3p_87(void); //                                                                                          
    void Form_I3p_88(void); //                                                                                          
    void Form_I3p_89(void); //                                                                                          
    void Form_I3p_90(void); //                                                                                          
    void Form_I3p_91(void); //                                                                                          
    void Form_I3p_92(void); //                                                                                          
    void Form_I3p_93(void); //                                                                                          
    void Form_I3p_94(void); //    
    void Form_I3p_95(void); //                                                                                          
    void Form_I3p_96(void); //                                                                                          
    void Form_I3p_97(void); //                                                                                          
    void Form_I3p_98(void); //                                                                                          
    void Form_I3p_99(void); //                                                                                          
    void Form_I3p_100(void); //                
                                                                                    
    void Form_I3p_101(void); //    
    void Form_I3p_102(void); // 
    void Form_I3p_103(void); //  
    void Form_I3p_104(void); //  
    void Form_I3p_105(void); //  
    void Form_I3p_106(void); //  
    void Form_I3p_107(void); //  
    void Form_I3p_108(void); //  
    void Form_I3p_109(void); //  
    void Form_I3p_110(void); //
    void Form_I3p_111(void); //                                   
    void Form_I3p_112(void); //                                   
    void Form_I3p_113(void); //                                   
    void Form_I3p_114(void); //                                   
    void Form_I3p_115(void); //                                   
    void Form_I3p_116(void); //                                   
    void Form_I3p_117(void); //                                   
    void Form_I3p_118(void); //   
    void Form_I3p_119(void); //                                   
    void Form_I3p_120(void); //                                   
    void Form_I3p_121(void); //                                   
    void Form_I3p_122(void); //                                   
    void Form_I3p_123(void); //                                   
    void Form_I3p_124(void); //                                   
    void Form_I3p_125(void); // 
    void Form_I3p_126(void); //     
    void Form_I3p_127(void); //                                                                         
    void Form_I3p_128(void); // 
    void Form_I3p_129(void); //                                                                         
    void Form_I3p_130(void); //                                                                         
    void Form_I3p_131(void); //                                                                         
    void Form_I3p_132(void); //                                                                         
    void Form_I3p_133(void); //                                                                         
    void Form_I3p_134(void); //    
    void Form_I3p_135(void); //                                                                         
    void Form_I3p_136(void); //                                                                         
    void Form_I3p_137(void); //                                                                         
    void Form_I3p_138(void); //                                                                         
    void Form_I3p_139(void); //                                                                         
    void Form_I3p_140(void); //                                                                         
    void Form_I3p_141(void); // 
    void Form_I3p_142(void); //                                                                                                                                                     
                                                                                                     
                                                                                                    
                                                                                
                                                                
                                                                       
                                                                                                                                                                                      
                                                                                                                      
                                                                                                                           
                                          


    void Form_I4e_1(void); // the fourth term first matrix
    void Form_I4p_2(void); //
    void Form_I4p_3(void); //
    void Form_I4p_4(void); //
    void Form_I4p_5(void); //
    void Form_I4p_6(void); //
    void Form_I4p_7(void); //
    void Form_I4p_8(void); //
    void Form_I4p_9(void); //
    void Form_I4p_10(void); //
    void Form_I4p_11(void); //
    void Form_I4p_12(void); //
    void Form_I4p_13(void); //
    void Form_I4p_14(void); //
    void Form_I4p_15(void); //
    void Form_I4p_16(void); //
    void Form_I4p_17(void); //
    void Form_I4p_18(void); //                            


    /* Functions for the second balance equation */
    void Form_IIJp_1(void);
    void Form_IIJp_2(void);
    void Form_IIJp_3(void);
    void Form_IIJp_4(void);
    void Form_IIJp_5(void);
    void Form_IIJp_6(void);
    void Form_IIJp_7(void);
    void Form_IIJp_8(void);
    /* Matrices from coupling*/
    void Form_IIJp_9(void);
    void Form_IIJp_10(void);
    void Form_IIJp_11(void);
    void Form_IIJp_12(void);
    void Form_IIJp_13(void);
    void Form_IIJp_14(void);
    void Form_IIJp_15(void);
    void Form_IIJp_16(void);                        


    void Form_II2e_1(void);
    void Form_II2p_2(void);
    void Form_II2p_3(void);
    void Form_II2p_4(void);
    void Form_II2p_5(void);
    void Form_II2p_6(void);
    void Form_II2p_7(void);
    void Form_II2p_8(void);
    void Form_II2p_9(void);
    /* Matrices from coupling*/    
    void Form_II2p_10(void);
    void Form_II2p_11(void);
    void Form_II2p_12(void);
    void Form_II2p_13(void);
    void Form_II2p_14(void);
    void Form_II2p_15(void);
    void Form_II2p_16(void);
    void Form_II2p_17(void);                            
        


    void Form_II3e_1(void);
    void Form_II3e_2(void);
    void Form_II3e_3(void);
    void Form_II3e_4(void);
    void Form_II3e_5(void);
    void Form_II3p_6(void);
    void Form_II3p_7(void);
    void Form_II3p_8(void);
    void Form_II3p_9(void);
    void Form_II3p_10(void);
    void Form_II3p_11(void);
    void Form_II3p_12(void);
    void Form_II3p_13(void);
    void Form_II3p_14(void);
    void Form_II3p_15(void);
    void Form_II3p_16(void);
    void Form_II3p_17(void);
    void Form_II3p_18(void);
    void Form_II3p_19(void);
    void Form_II3p_20(void);
    void Form_II3p_21(void);
    void Form_II3p_22(void);
    void Form_II3p_23(void);
    void Form_II3p_24(void);
    void Form_II3p_25(void);
    void Form_II3p_26(void);
    void Form_II3p_27(void);
    void Form_II3p_28(void);
    void Form_II3p_29(void);
    void Form_II3p_30(void);
    void Form_II3e_31(void);
    void Form_II3e_32(void);
    void Form_II3e_33(void);
    void Form_II3p_34(void);
    void Form_II3p_35(void);
    void Form_II3p_36(void);
    void Form_II3p_37(void);
    void Form_II3p_38(void);
    void Form_II3p_39(void);
    void Form_II3p_40(void);
    void Form_II3p_41(void);
    void Form_II3p_42(void);
    void Form_II3p_43(void);
    void Form_II3p_44(void);
    void Form_II3p_45(void);
    void Form_II3p_46(void);
    void Form_II3p_47(void);
    void Form_II3p_48(void);
    /* Matrices from Del(delgammachi) */
    void Form_II3p_49(void);  
    void Form_II3p_50(void);    
    void Form_II3p_51(void);    
    void Form_II3p_52(void);    
    void Form_II3p_53(void);    
    void Form_II3p_54(void);    
    void Form_II3p_55(void);
    void Form_II3p_56(void);                                                        
    void Form_II3p_57(void);  
    void Form_II3p_58(void);    
    void Form_II3p_59(void);    
    void Form_II3p_60(void);    
    void Form_II3p_61(void);    
    void Form_II3p_62(void);    
    void Form_II3p_63(void);
    void Form_II3p_64(void);    
    void Form_II3p_65(void);  
    void Form_II3p_66(void);    
    void Form_II3p_67(void);    
    void Form_II3p_68(void);    
    void Form_II3p_69(void);    
    void Form_II3p_70(void);    
    void Form_II3p_71(void);
    void Form_II3p_72(void);    
    /* Matrices from coupling */
    void Form_II3p_73(void);    
    void Form_II3p_74(void);    
    void Form_II3p_75(void);    
    void Form_II3p_76(void);    
    void Form_II3p_77(void);    
    void Form_II3p_78(void);    
    void Form_II3p_79(void);    
    void Form_II3p_80(void);
    void Form_II3p_81(void);    
    void Form_II3p_82(void);    
    void Form_II3p_83(void);    
    void Form_II3p_84(void);    
    void Form_II3p_85(void);    
    void Form_II3p_86(void);    
    void Form_II3p_87(void);    
    void Form_II3p_88(void);                                        
                                    




    void Form_II4e_1(void);
    void Form_II4p_2(void);
    void Form_II4p_3(void);
    void Form_II4p_4(void);
    void Form_II4p_5(void);
    void Form_II4p_6(void);
    void Form_II4p_7(void);
    void Form_II4p_8(void);
    void Form_II4p_9(void);

    /////////////////////////////////////////////////////////



    //////////////////////////////////////////////////////////
    /////FUNCTIONS FINISH HERE FOR MICROMORPHIC MATRICES////
    //////////////////////////////////////////////////////////

protected:

    /** extract natural boundary condition information */
    void TakeNaturalBC(const ParameterListT& list);

    /** apply traction boundary conditions to displacement equations */
    void ApplyTractionBC(void);

    /** update traction BC data for displacement equations */
    void SetTractionBC(void);

    /* traction data */
    ArrayT<Traction_CardT> fTractionList;
    int fTractionBCSet;

    /** \name arrays with local ordering */
    /*@{*/
    LocalArrayT fLocInitCoords;   /**< initial coords with local ordering */
    LocalArrayT fLocDisp;              /**< solid displacements with local ordering  */
    /*@}*/

    /** \name work space */
    /*@{*/
    dArrayT fNEEvec; /**< work space vector: [element DOF] */
    dArrayT fDOFvec; /**< work space vector: [nodal DOF]   */
    /*@}*/

};


inline const ShapeFunctionT& FSMicromorphic3DT::ShapeFunctionDispl(void) const
{
#if __option(extended_errorcheck)
    if (!fShapes_displ)
        ExceptionT::GeneralFail("FSMicromorphic3DT::ShapeFunctionDispl", "no displ shape functions");
#endif
    return *fShapes_displ;
}

inline const ShapeFunctionT& FSMicromorphic3DT::ShapeFunctionMicro(void) const
{
#if __option(extended_errorcheck)
    if (!fShapes_micro)
        ExceptionT::GeneralFail("FSMicromorphic3DT::ShapeFunctionMicro", "no micro shape functions");
#endif
    return *fShapes_micro;
}

/* return the geometry code */
inline GeometryT::CodeT FSMicromorphic3DT::GeometryCode(void) const
{ return fGeometryCode_displ; }


} // namespace Tahoe
#endif /* _FS_MICROMORPHIC_3D_T_H_ */



