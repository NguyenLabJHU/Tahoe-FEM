/*
  File GradientTools.cpp
*/

#include "GradientTools.h"
#include "Utils.h"

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

GradientTools::GradientTools(int numip, int numnodes, int numsd) :
  fNumIP       (numip),
  fNumNodes    (numnodes),  // only vertex nodes: 4 (Quad), 8 (Hexa)
  fNumSD       (numsd),
  fNa          (fNumIP, fNumNodes),
  fLDNa        (fNumIP),
  fGDNa        (fNumIP),
  fNodalExtrap (fNumNodes, fNumIP),
  fWeights     (fNumIP),
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
    if (fNumNodes != 4)
       throwRunTimeError("GradientTools::Constructor: numnode != 4 (nsd=2)");
    if (fNumIP != 1 &&
        fNumIP != 4 && 
        fNumIP != 9)
       throwRunTimeError("GradientTools::Constructor: numint != 1, 4 or 9 (nsd=2)");
  }
  else {    // (fNumSD == 3)
    if (fNumNodes != 8)
      throwRunTimeError("GradientTools::Constructor: numnodes != 8 (nsd=3)");
    if (fNumIP != 1 &&
        fNumIP != 8 &&
        fNumIP != 27)
      throwRunTimeError("GradientTools::Constructor: numint != 1, 8 or 27 (nsd=3)");
  }

  // addional allocation for shape function derivatives
  for (int i = 0; i < fNumIP; i++) {
     fLDNa[i].Dimension(fNumSD, fNumNodes);
     fGDNa[i].Dimension(fNumSD, fNumNodes);
  }

  // additional allocation for spatial gradient of a tensor T
  for (int i = 0; i < kNSD; i++)
     fGradT[i].Dimension(kNSD, kNSD); 

  //initialize arrays
  fNa = 0.0;
  fNodalExtrap = 0.0;
  for (int i = 0; i < fNumIP; i++) fLDNa[i] = 0.0;

  // set local shape functions and their derivatives
  if (fNumSD == 2) SetLocalShape2D();     // only vertex nodes: 4 (Quad)
  if (fNumSD == 3) SetLocalShape3D();     // only vertex nodes: 8 (Hexa)

  // set extrapolation matrix of ipvariables to npvariables
  SetNodalExtrapolation();
}

// dNa/dE (LDNa) -> dNa/dX (GDNa): coords = Initial_X
// dNa/dE (LDNa) -> dNa/dx (GDNa): coords = Current_X
/* coords (#nodes x #sd), fLDNa[#ip x (#sd x #nodes)], fGDNa[#ip x (#sd x #nodes)] */
/* note: #ip=fLDNa.Length(), #sd=fLDNa.MajorDim() */
void GradientTools::ComputeGDNa(const LocalArrayT& coords)
{
  // loop over all integration points
  for (int i = 0; i < fNumIP; i++)
    {
      // compute Jacobian matrix
      Jacobian(coords, fLDNa[i]); 
      double det = fJac.Det();
      
      // element check
      if (det < kSmall) throwRunTimeError("GradientTools::ComputeGDNa: det(Jac) < ksmall");

      dMatrixT& jac_inv = fJac.Inverse();

      // global shape function derivatives
      ChangeOfVariables(jac_inv, fLDNa[i], fGDNa[i]);
    }
}

/* npvalues(#nodes x (#sd x #sd)), ipvalues(#IP x (#sd x #sd)), Na(#IP x #nodes) */
/* note: #ip=Na.MajorDim(), #nodes=Na.MinorDim() */
void GradientTools::InterpolateTensor(const ArrayT<dMatrixT>& npvalues, 
                                      ArrayT<dMatrixT>& ipvalues)
{
  // interpolate to all integration points
  for (int i = 0; i < fNumIP; i++)
    {
      ipvalues[i] = 0.0;
      for (int j = 0; j < fNumNodes; j++)
	ipvalues[i].AddScaled(fNa(i, j), npvalues[j]);
    }
}

/* ipvalues[#IP x (#sd x #sd)], npvalues[#nodes x (#sd x #sd)], Ma(#nodes x #IP) */
/* note: #nodes=Ma.MajorDim(), #ip=Ma.MinorDim() */
void GradientTools::ExtrapolateTensor(const ArrayT<dMatrixT>& ipvalues, 
                                      ArrayT<dMatrixT>& npvalues)
{
  // extrapolate to all nodal points
  for (int i = 0; i < fNumNodes; i++)
    {
      npvalues[i] = 0.0;
      for (int j = 0; j < fNumIP; j++)
	npvalues[i].AddScaled(fNodalExtrap(i, j), ipvalues[j]);
    }
}

// spatial gradient of a tensor (dA/dx) at particular IP
/* nodal [#nodes x (#sd x #sd)], tensorGrad [#knsd x (#sd x #sd)] */
/* fGDNa[#IP x (#sd x #nodes)] */
void GradientTools::GradientOfTensorAtIP(const ArrayT<dMatrixT>& nodal, 
   ArrayT<dMatrixT>& tensorGrad, int ip)
{
  //dimensions of tensor grad (knsd = 3; a must!)
  int knsd = tensorGrad.Length();

  // check dimensions (needed when tensorGrad is external variable)
  if (knsd != kNSD)
     throwRunTimeError("GradientTools::GradientOfTensorAtIP: knsd != kNSD");

  // intialize spatial gradient d()/dx
  for (int i = 0; i < knsd; i++) tensorGrad[i] = 0.;

  // compute dFe/dx at integration point
  if (fNumSD == 2)
    {
      double* dx1 = fGDNa[ip](0);
      double* dx2 = fGDNa[ip](1);

      for (int i = 0; i < fNumNodes; i++)
        {
          tensorGrad[0].AddScaled(*dx1, nodal[i]);
          tensorGrad[1].AddScaled(*dx2, nodal[i]);

          dx1++; dx2++;
        }
    }
  else // (fNumSD ==3)
    {
      double* dx1 = fGDNa[ip](0);
      double* dx2 = fGDNa[ip](1);
      double* dx3 = fGDNa[ip](2);

      for (int i = 0; i < fNumNodes; i++)
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
/* Tnodal [#nodes x (knsd x knsd)],  CurlT (knsd x knsd) */
/* note: for X^p -> Tnodal=Fp;;; for X^e -> Tnodal=Fe^(-1) */
void GradientTools::CurlOfTensorAtIP(const ArrayT<dMatrixT>& Tnodal, 
   dMatrixT& CurlT, int ip)
{
  //dimensions of tensor grad
  int knsd = CurlT.Rows();

  // check dimensions (fGradT is local variable)
  if (knsd != fGradT.Length())
     throwRunTimeError("GradientTools::CurlOfTensorAtIP: knsd != fGraT.Length()");

  // gradient of tensor T
  GradientOfTensorAtIP(Tnodal, fGradT, ip);

  // compute CurlT
  for (int Jbar = 0; Jbar < knsd; Jbar++)
     {
        CurlT(0,Jbar) = fGradT[1](Jbar,2) - fGradT[2](Jbar,1);
        CurlT(1,Jbar) =-fGradT[0](Jbar,2) + fGradT[2](Jbar,0);
        CurlT(2,Jbar) = fGradT[0](Jbar,1) - fGradT[1](Jbar,0);
     }
}

/* PRIVATE MEMBER FUNCTIONS */

// set shape functions and their derivatives at IPs for 2D
/* fNa(#IP x #nodes), fLDNa[#IP x (#sd x #nodes)] */
void GradientTools::SetLocalShape2D()
{
    // integration point coordinates
    double  ra[9] = {-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0,-1.0, 0.0};
    double  sa[9] = {-1.0,-1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 0.0, 0.0};
    double *xa, *ya;
    double g;
  
    // integration weights
    switch (fNumIP)
    {
      case 1:	
	xa = ra;			  
	ya = sa;			  
	g = 0.0;
	fWeights[0] = 4.0;
	break;
	
      case 4:
	xa = ra;			  
	ya = sa;			  
	g = 1.0/sqrt3;
	fWeights[0] = 1.0;
	fWeights[1] = 1.0;
	fWeights[2] = 1.0;
	fWeights[3] = 1.0;
	break;
	
      case 9:
        {
	xa = ra;			  
	ya = sa;			  
	g = sqrt(3.0/5.0);
	double a = 25.0/81.0;
	double b = 40.0/81.0;
	double c = 64.0/81.0;
	fWeights[0] = a;
	fWeights[1] = a;
	fWeights[2] = a;
	fWeights[3] = a;
	fWeights[4] = b;
	fWeights[5] = b;
	fWeights[6] = b;
	fWeights[7] = b;
	fWeights[8] = c;
	break;
        }
	
      default:
	throwRunTimeError("Utils::SetLocalShape: Bad numint (nsd=2)");
    }
  
    // shape functions and derivatives at IPs (using vertex nodes)
    for (int in = 0; in < fNumIP; in++)	
      {	
	double* na  = fNa(in);
	double* nax = fLDNa[in](0);
	double* nay = fLDNa[in](1);
	
	double r = g*xa[in];
	double s = g*ya[in];
	
	for (int lnd = 0; lnd < fNumNodes; lnd++)
	  {
	    double tempr1 = 1.0 + ra[lnd]*r;
	    double temps1 = 1.0 + sa[lnd]*s;
	    
	    *na++  = 0.25*tempr1*temps1;
	    *nax++ = 0.25*ra[lnd]*temps1;
	    *nay++ = 0.25*tempr1*sa[lnd];
	  }
      }
}

// set shape functions and their derivatives at IPs for 3D
/* fNa(#IP x #nodes), fLDNa[#IP x (#sd x #nodes)] */
void GradientTools::SetLocalShape3D()
{
    /* integration point coordinates */
    double  ra[] = {-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0};
    double  sa[] = {-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0};
    double  ta[] = {-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0};
    double  xa_27[27], ya_27[27], za_27[27];
    double  *xa, *ya, *za;
    double  g;
    
    /* integration weights */
    switch (fNumIP)
      {
      case 1:	
	g = 0.0;
	xa = ra;
	ya = sa;
	za = ta;
	fWeights[0] = 8.0;
	break;
	
      case 8:
	g = 1.0/sqrt3;
	xa = ra;
	ya = sa;
	za = ta;
	fWeights = 1.0;
	break;
	
      case 27:
        {
	/* coordinates */
	double b1 = sqrt(3.0/5.0);
	double b_1D[3] = {-b1, 0.0, b1}; 
	
	/* weights */
	double w1 = 5.0/9.0;
	double w2 = 8.0/9.0;
	double w_1D[3] = {w1, w2, w1};
	int x_i = 0;
	int y_i = 0;
	int z_i = 0;
	for (int i = 0; i < 27; i++)
	  {
	    xa_27[i]   = b_1D[x_i];
	    ya_27[i]   = b_1D[y_i];
	    za_27[i]   = b_1D[z_i];
	    fWeights[i] = w_1D[x_i]*w_1D[y_i]*w_1D[z_i];
	    
	    if (++x_i == 3)
	      {
		x_i = 0;
		if (++y_i == 3)
		  {
		    y_i = 0;
		    z_i++;
		  }
	      }
	  }						
	
	xa = xa_27;
	ya = ya_27;
	za = za_27;
	g  = 1.0;		
	break;
        }
	
      default:
	throwRunTimeError("Utils::SetLocalShape: Bad numint (nsd=3)");
      }
    
    // shape functions and derivatives at IPs (using vertex nodes)
    for (int in = 0; in < fNumIP; in++)	
      {
	double* na  = fNa(in);
	double* nax = fLDNa[in](0);
	double* nay = fLDNa[in](1);
	double* naz = fLDNa[in](2);
	
	double r = g*xa[in];
	double s = g*ya[in];
	double t = g*za[in];
	
	for (int lnd = 0; lnd < fNumNodes; lnd++)
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
}

// set extrapolation matrix of ipvariables to npvariables 
/* fNodalExtrap(#nodes x #IP), fNa(#IP x #nodes) */
void GradientTools::SetNodalExtrapolation()
{
  // nodal extrapolation matrix (least square smoothing)
  dArrayT  v1 (fNumNodes);
  dArrayT  v2 (fNumNodes);
  dMatrixT NaxNaT     (fNumNodes);
  dMatrixT NaxNaT_inv (fNumNodes);

  NaxNaT = 0.0;
  for (int i = 0; i < fNumIP; i++)
    {
      fNa.RowCopy(i, v1);
      NaxNaT_inv.Outer(v1, v1);    // NaxNaT_inv is dummy here
      NaxNaT.AddScaled(fWeights[i], NaxNaT_inv);
    }

  NaxNaT_inv = MatrixInversion(NaxNaT);

  for (int i = 0; i < fNumIP; i++)
    {
      fNa.RowCopy(i, v1);
      NaxNaT_inv.Multx(v1, v2);
      v2 *= fWeights[i];
      fNodalExtrap.SetColumn(i, v2);
    }
}

/* nodal (#nodes x #sd), LDNa (#sd x #nodes), fJac (#sd x #sd) */
/* note: #sd=LDNa.MajorDim(), #nodes=LDNa.MinorDim() */
void GradientTools::Jacobian(const LocalArrayT& nodal, const dArray2DT& LDNa)
{
  // jacobian
  fJac = 0.;
  for (int i = 0; i < fNumNodes; i++)
    {
      LDNa.ColumnCopy(i, fArray2);
      for (int j = 0; j < fNumSD; j++) fArray1[j] = nodal(i,j);
      fMatx1.Outer(fArray1, fArray2);
      fJac += fMatx1;
    }
}

// each IP at a time
/* dNa/dE -> dNa/dX : tensor = Jac^(-1) */
/* dNa/dX -> dNa/dx : tensor = F^(-1)   */
/* tensor (#sd x #sd), DNa (#sd x #nodes) */
/* note: #sd=DNa.MajorDim(), #nodes=DNa.MinorDim() */
void GradientTools::ChangeOfVariables(const dMatrixT& tensor, const dArray2DT& DNa_known, 
                                      dArray2DT& DNa_unknown)
{
  // change of variables in derivatives
  for (int i = 0; i < fNumNodes; i++)
    {
      DNa_known.ColumnCopy(i, fArray1);
      tensor.MultTx(fArray1, fArray2);
      DNa_unknown.SetColumn(i, fArray2);
    }
}
