/*
  File GradientTools_C.cpp
*/

#include "GradientTools_C.h"
#include "Utils.h"
#include "ifstreamT.h"
#include "StringT.h"
#include "ExceptionT.h"
#include "toolboxConstants.h"

#include "LocalArrayT.h"
#include "dArrayT.h"
#include "iArrayT.h"
#include "dArray2DT.h"
#include "dMatrixT.h"


using namespace Tahoe;

const int     kNSD = 3;
const double sqrt3 = sqrt(3.0);

GradientTools_C::GradientTools_C(int numip, int numsd) :
  fNumIP       (numip),
  fNumSD       (numsd),
  fNa          (1,fNumIP),
  fArray1      (fNumSD),
  fArray2      (fNumSD),
  fJac         (fNumSD),
  fMatx1       (fNumSD),
  fGradT       (kNSD)
{
  // check spatial dimensions
  if (fNumSD != 2 && fNumSD != 3)
     throwRunTimeError("GradientTools::Constructor: nsd != 2 or 3");

  // additional dimension checks
  if (fNumSD == 2) {
    if (fNumIP != 4)
       throwRunTimeError("GradientTools::Constructor: numint != 4 (nsd=2)");
  }
  else {    // (fNumSD == 3)
    if (fNumIP != 8)
      throwRunTimeError("GradientTools::Constructor: numint != 8 (nsd=3)");
  }

  // addional allocation for shape function derivatives
  fLDNa.Dimension(fNumSD, fNumIP);
  fGDNa.Dimension(fNumSD, fNumIP);

  // additional allocation for spatial gradient of a tensor T
  for (int i = 0; i < kNSD; i++)
     fGradT[i].Dimension(kNSD, kNSD); 

  //initialize arrays
  fNa = 0.0;
  fLDNa = 0.0;

  // set local shape functions and their derivatives
  if (fNumSD == 2) SetLocalShape2D();     // 4 IPs (Quad)
  if (fNumSD == 3) SetLocalShape3D();     // 8 IPs (Hexa)
}

// dNa/dE (LDNa) -> dNa/dX (GDNa): coords = Initial_X
// dNa/dE (LDNa) -> dNa/dx (GDNa): coords = Current_X
/* coords (#ip x #sd), fLDNa(#sd x #ip), fGDNa(#sd x #ip) */
/* note: #ip=fLDNa.MinorDim(), #sd=fLDNa.MajorDim() */
void GradientTools_C::ComputeGDNa(const LocalArrayT& coords)
{
  // compute Jacobian matrix
  Jacobian(coords, fLDNa); 

  // element check
  double det = fJac.Det();
  if (det < kSmall) throwRunTimeError("GradientTools::ComputeGDNa: det(Jac) < ksmall");

  dMatrixT& jac_inv = fJac.Inverse();

  // global shape function derivatives
  ChangeOfVariables(jac_inv, fLDNa, fGDNa);
}

/* ipvalues(#ip x (#sd x #sd)), cvalue(#sd x #sd), Na(1 x #IP) */
/* note: #ip=Na.MinorDim() */
void GradientTools_C::InterpolateTensorAtCenter(const ArrayT<dMatrixT>& ipvalues, 
                                                dMatrixT& cvalue)
{
  // interpolate to center of element
  cvalue = 0.0;
  for (int j = 0; j < fNumIP; j++)
    cvalue.AddScaled(fNa(0,j), ipvalues[j]);
}

// spatial gradient of a tensor (dA/dx) at particular IP
/* nodal [#ip x (#sd x #sd)], tensorGrad [#knsd x (#sd x #sd)] */
/* fGDNa(#sd x #ip) */
void GradientTools_C::GradientOfTensorAtCenter(const ArrayT<dMatrixT>& nodal, 
   ArrayT<dMatrixT>& tensorGrad)
{
  //dimensions of tensor grad (knsd = 3; a must!)
  int knsd = tensorGrad.Length();

  // check dimensions (needed when tensorGrad is external variable)
  if (knsd != kNSD)
     throwRunTimeError("GradientTools::GradientOfTensorAtIP: knsd != kNSD");

  // intialize spatial gradient d()/dx
  for (int i = 0; i < knsd; i++) tensorGrad[i] = 0.;

  // compute dFe/dx at center of element
  if (fNumSD == 2)
    {
      double* dx1 = fGDNa(0);
      double* dx2 = fGDNa(1);

      for (int i = 0; i < fNumIP; i++)
        {
          tensorGrad[0].AddScaled(*dx1, nodal[i]);
          tensorGrad[1].AddScaled(*dx2, nodal[i]);

          dx1++; dx2++;
        }
    }
  else // (fNumSD ==3)
    {
      double* dx1 = fGDNa(0);
      double* dx2 = fGDNa(1);
      double* dx3 = fGDNa(2);

      for (int i = 0; i < fNumIP; i++)
        {
          tensorGrad[0].AddScaled(*dx1, nodal[i]);
          tensorGrad[1].AddScaled(*dx2, nodal[i]);
          tensorGrad[2].AddScaled(*dx3, nodal[i]);

          dx1++; dx2++; dx3++;
        }
    }
}

// Compute the Curl (or curl) of a Tensor: A^T
// (important: see order of subscripts)
// 1. Fundamental Plastic Dislocation Tensor
//          X^p =  Curl(A)^T = Nabla_X x A^T
//          A maps: X->Xbar ; i.e. A=Fp
//    then:
//        Curl(A)^T = d/dX_I e_I x (A_JbarJ e_Jbar e_J)^T
//                  = E_IJK dA_JbarJ/dX_I e_K e_Jbar
//        => [Curl(A)^T]_KJbar = E_IJK dA_JbarJ/dX_I 
//    (applicable also for maps: X->x; then Jbar=j)
//
// 2. True Dislocation Density Tensor
//          X^e = curl(A)^(-T) = Nabla_x x A^(-T)
//          A^(-1) maps: x->Xbar ; e.g. A=Fe
//    then: 
//        curl(A)^(-T) = d/dx_i e_i x (A_Jbarj^(-1) e_Jbar e_j)^T
//                     = E_ijk dA_Jbarj^(-1)/dx_i e_k e_Jbar
//          => [Curl(A)^(-T)]_kJbar = E_ijk dA_Jbarj^(-1)/dx_i
//    (applicable also for maps: x->X; then Jbar=J)
//
/* Tnodal [#ip x (knsd x knsd)],  CurlT (knsd x knsd) */
/* note: for X^p -> Tnodal=Fp;;; for X^e -> Tnodal=Fe^(-1) */
void GradientTools_C::CurlOfTensorAtCenter(const ArrayT<dMatrixT>& Tnodal, 
   dMatrixT& CurlT)
{
  //dimensions of tensor grad
  int knsd = CurlT.Rows();

  // check dimensions (fGradT is local variable)
  if (knsd != fGradT.Length())
     throwRunTimeError("GradientTools::CurlOfTensorAtIP: knsd != fGraT.Length()");

  // gradient of tensor T
  GradientOfTensorAtCenter(Tnodal, fGradT);

  // compute CurlT
  for (int Jbar = 0; Jbar < knsd; Jbar++)
     {
        CurlT(0,Jbar) = fGradT[1](Jbar,2) - fGradT[2](Jbar,1);
        CurlT(1,Jbar) =-fGradT[0](Jbar,2) + fGradT[2](Jbar,0);
        CurlT(2,Jbar) = fGradT[0](Jbar,1) - fGradT[1](Jbar,0);
     }
}

/* PRIVATE MEMBER FUNCTIONS */

// set shape functions and their derivatives at center of IP element - 2D
/* fNa(1 x #IP), fLDNa(#sd x #IP); #IP=4 */
void GradientTools_C::SetLocalShape2D()
{
    // nodal coordinates of IP element
    double  ra[4] = {-1.0, 1.0, 1.0,-1.0};
    double  sa[4] = {-1.0,-1.0, 1.0, 1.0};
  
    // coordinates of center of IP element
    double r = 0.;
    double s = 0.;
	
    // shape functions and derivatives at center of IP element
    double* na  = fNa(0);
    double* nax = fLDNa(0);
    double* nay = fLDNa(1);
	
    for (int lnd = 0; lnd < fNumIP; lnd++)
       {
         double tempr1 = 1.0 + ra[lnd]*r;
         double temps1 = 1.0 + sa[lnd]*s;
    
         *na++  = 0.25*tempr1*temps1;
         *nax++ = 0.25*ra[lnd]*temps1;
         *nay++ = 0.25*tempr1*sa[lnd];
       }
}

// set shape functions and their derivatives at center of IP element - 3D
/* fNa(1 x #IP), fLDNa(#sd x #nodes); #IP=8 */
void GradientTools_C::SetLocalShape3D()
{
    /* integration point coordinates */
    double  ra[8] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
    double  sa[8] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
    double  ta[8] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};
    
    // coordinates of center of element
    double r = 0.;
    double s = 0.;
    double t = 0.;

    // shape functions and derivatives at center of element
    double* na  = fNa(0);
    double* nax = fLDNa(0);
    double* nay = fLDNa(1);
    double* naz = fLDNa(2);
	
    for (int lnd = 0; lnd < fNumIP; lnd++)
       {
         double tempr1 = 1.0 + ra[lnd]*r;
         double temps1 = 1.0 + sa[lnd]*s;
         double tempt1 = 1.0 + ta[lnd]*t;
    
         *na++  = 0.125*tempr1*temps1*tempt1;
         *nax++ = 0.125*ra[lnd]*temps1*tempt1;
         *nay++ = 0.125*tempr1*sa[lnd]*tempt1;
         *naz++ = 0.125*tempr1*temps1*ta[lnd];
       }
}

/* nodal (#ip x #sd), LDNa (#sd x #ip), fJac (#sd x #sd) */
/* note: #sd=LDNa.MajorDim(), #ip=LDNa.MinorDim() */
void GradientTools_C::Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNa)
{
  // jacobian
  fJac = 0.;
  for (int i = 0; i < fNumIP; i++)
    {
      LDNa.ColumnCopy(i, fArray2);
      for (int j = 0; j < fNumSD; j++) fArray1[j] = nodal(i,j);
      fMatx1.Outer(fArray1, fArray2);
      fJac += fMatx1;
    }
}

/* dNa/dE -> dNa/dX : tensor = Jac^(-1) */
/* dNa/dX -> dNa/dx : tensor = F^(-1)   */
/* tensor (#sd x #sd), DNa (#sd x #ip) */
/* note: #sd=DNa.MajorDim(), #ip=DNa.MinorDim() */
void GradientTools_C::ChangeOfVariables(const dMatrixT& tensor, 
                                        const dArray2DT& DNa_known, 
                                        dArray2DT& DNa_unknown)
{
  // change of variables in derivatives
  for (int i = 0; i < fNumIP; i++)
    {
      DNa_known.ColumnCopy(i, fArray1);
      tensor.MultTx(fArray1, fArray2);
      DNa_unknown.SetColumn(i, fArray2);
    }
}
