/* $Id: SmallStrainQ2P1.cpp,v 1.6 2004-07-15 08:28:31 paklein Exp $ */
#include "SmallStrainQ2P1.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "MaterialListT.h"
#if 0
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"
#endif
#include "SSMatSupportT.h"

using namespace Tahoe;

/* constructor */
SmallStrainQ2P1::SmallStrainQ2P1(const ElementSupportT& support, const FieldT& field):
  SmallStrainT(support),
  fpdof(NumSD()+1),
  fthird(1.0/3.0)
{}

/* called immediately after constructor */
void SmallStrainQ2P1::Initialize(void)
{
ExceptionT::GeneralFail("SmallStrainQ2P1::Initialize", "out of date");
#if 0	
  /* inherited */
  SolidElementT::Initialize();
  
  int nsd = (NumSD()==2) ?4:NumSD();
  int nel = NumElements();
  int nen = NumElementNodes();
  int nip = NumIP();
  
  /* check geometry code and number of element nodes -> Q1 */
  if (GeometryCode() == GeometryT::kQuadrilateral) {
    if (nen < 9 ) {
      cout << "\n SmallStrainQ2P1::Initialize: expecting 9 node quad: " 
	   << NumElementNodes() << endl;
      throw ExceptionT::kBadInputValue;
    }   
  }
  else {
    cout << "\n SmallStrainQ2P1::Initialize: expecting quad geometry:"
	 << GeometryCode() << endl;
    throw ExceptionT::kBadInputValue;
  }	
  if (nsd ==4)
    fD.Dimension(nsd);

  /*initializes Bbar and utility tensors*/	
  fBbar.Dimension(dSymMatrixT::NumValues(nsd), NumSD()*nen);
  fB_dev.Dimension(dSymMatrixT::NumValues(nsd), NumSD()*nen);
  fBbar_dil.Dimension(dSymMatrixT::NumValues(nsd), NumSD()*nen);
  fGradbar.Dimension(nip);
  for (int i=0; i<nip; i++)
    fGradbar[i].Dimension(NumSD(),nen);
  fGradTranspose.Dimension(NumSD());
  /* what's needed */
  bool need_strain = false;
  bool need_strain_last = false;
  for (int i = 0; i < fMaterialNeeds.Length(); i++){
    const ArrayT<bool>& needs = fMaterialNeeds[i];
    need_strain = need_strain || needs[fNeedsOffset + kstrain];
    need_strain_last = need_strain_last || needs[fNeedsOffset + kstrain_last];
  }

  /* allocate deformation gradient list */
  if (need_strain) {
    fStrain_List.Dimension(NumIP());
    for (int i = 0; i < NumIP(); i++)
      fStrain_List[i].Dimension(nsd);
    fTheta_List.Dimension(nel,nip);
    fTheta_List = 0.0;
  }
	
  /* allocate "last" deformation gradient list */
  if (need_strain_last) {
    fStrain_last_List.Dimension(NumIP());
    for (int i = 0; i < NumIP(); i++)
      fStrain_last_List[i].Dimension(nsd);
    fTheta_List_last.Dimension(nel,nip);
    fTheta_List_last = 0.0;
  }
  
  fPressure_List.Dimension(nel,nip);
  fPressure_List = 0.0;
  
  fMShapes.Dimension(nip,fpdof);
  fGamma.Dimension(fpdof);
  fip_coords.Dimension(NumSD());
  fH_inv.Dimension(fpdof);
  fAVec.Dimension(fpdof);
  fAMat.Dimension(fpdof,NumSD());
  fBMat.Dimension(fpdof,NumSD());
#endif
}

void SmallStrainQ2P1::CloseStep(void)
{
  /*inherited*/
  ContinuumElementT::CloseStep();
  fTheta_List_last = fTheta_List;
}

GlobalT::RelaxCodeT SmallStrainQ2P1::ResetStep(void)
{
	fTheta_List = fTheta_List_last;

	/* inherited */
	return ContinuumElementT::ResetStep();
}

void SmallStrainQ2P1::ReadRestart(istream& in)
{
  /*inherited*/
  ContinuumElementT::ReadRestart(in);
  in >> fTheta_List;
  fTheta_List_last = fTheta_List;
}

void SmallStrainQ2P1::WriteRestart(ostream& out) const
{
  /*inherited*/
  ContinuumElementT::WriteRestart(out);
  out << fTheta_List; 
}


/***********************************************************************
* Protected
***********************************************************************/
void SmallStrainQ2P1::SetLocalArrays(void)
{
  /*inherited*/
  SmallStrainT::SetLocalArrays();
  fLocDispTranspose.Dimension(fLocDisp.Length());
}
/* compute the measures of strain/deformation over the element */
void SmallStrainQ2P1::SetGlobalShape(void)
{
  /* inherited */
  SolidElementT::SetGlobalShape();
  
  /* material information */
  int matnum = CurrentElement().MaterialNumber();
  const ArrayT<bool>& needs = fMaterialNeeds[matnum];
  
  /* element information */
  int nen = NumElementNodes();
  int nsd = (NumSD()==2)? 4: NumSD();
  int nip = fShapes->NumIP();
  int elem = CurrElementNumber();
  
  const double* jac = fShapes->IPDets();
  const double* W = fShapes->IPWeights();  
  fH_inv = 0.0;
  
  double *pshapes = fMShapes.Pointer();
  fShapes->TopIP();
  while (fShapes->NextIP()) {
    /* pressure shape functions */
    *pshapes++ = 1.0;
    fGamma[0] = 1.0;
    IP_Coords(fip_coords);
    for (int i = 0; i<fip_coords.Length(); i++) {
      *pshapes++ = fip_coords[i];
      fGamma[i+1] = fip_coords[i];
    }
    /* mass matrix formed from pressure shape functions */
    fH_inv.Outer(fGamma,fGamma, (*W++)*(*jac++), dMatrixT::kAccumulate);
  }
  fH_inv.Inverse();
  /*set shape function gradients for mixed dilation*/
  Set_Gradbar();
  
  /* set mixed strains */
  double* ptheta = fTheta_List(elem);
  double* ptheta_last = fTheta_List_last(elem);
  for(int ip = 0; ip < nip; ip++){	
      Set_Bbar(fShapes->Derivatives_U(ip), fGradbar[ip], fBbar, fB_dev, fBbar_dil);
      /* deformation gradient */
      if (needs[fNeedsOffset + kstrain]) {
	/* transpose displacement array */
	fLocDisp.ReturnTranspose(fLocDispTranspose);

	/*compute strain using Bbar*/
	dSymMatrixT& strain = fStrain_List[ip];
	fBbar.Multx(fLocDispTranspose, strain);
	strain.ScaleOffDiagonal(0.5);
	*ptheta++ = strain.Trace();
      }

      /* "last" deformation gradient */
      if (needs[fNeedsOffset + kstrain_last]) {
	  /* transpose displacement array */
	  fLocLastDisp.ReturnTranspose(fLocDispTranspose);
		 
	  /*compute strain using Bbar*/
	  dSymMatrixT& strain = fStrain_last_List[ip];
	  fBbar.Multx(fLocDispTranspose, strain);
	  strain.ScaleOffDiagonal(0.5);
	  *ptheta_last++ = strain.Trace();
      }
  }	      
}

void SmallStrainQ2P1::Set_Gradbar(void)
{
  int nsd = (NumSD()==2)?4:3;
  for (int a = 0; a<NumElementNodes(); a++) {
    const double* jac = fShapes->IPDets();
    const double* W = fShapes->IPWeights();  
    fAMat = 0.0;
    fBMat = 0.0;

    fShapes->TopIP();
    while (fShapes->NextIP()) {
      double scale = (*jac++)*(*W++); 
      const dArray2DT& DNa = fShapes->Derivatives_U();	       
      fMShapes.RowCopy(CurrIP(),fGamma);
      DNa.ColumnCopy(a,fGradTranspose);
      fAMat.Outer(fGamma, fGradTranspose, scale, dMatrixT::kAccumulate);
    }
    fBMat.MultAB(fH_inv, fAMat, dMatrixT::kWhole);
     for (int i = 0; i<NumIP(); i++)  {
       fMShapes.RowCopy(i,fGamma);
       fBMat.MultTx(fGamma, fGradTranspose);
       fGradbar[i].SetColumn(a, fGradTranspose);
     }
  }
}

void SmallStrainQ2P1::Set_Bbar(const dArray2DT& DNa, const dArray2DT& DNabar, 
			       dMatrixT& Bbar,dMatrixT& B_dev, dMatrixT& Bbar_dil)
{
  double* pBdev = B_dev.Pointer();
  double* pBdil = Bbar_dil.Pointer();

  if (DNabar.MajorDim()==2){
    const double* pDx = DNa(0); const double* pDbx = DNabar(0); 
    const double* pDy = DNa(1); const double* pDby = DNabar(1);  
  
   for (int a = 0; a<NumElementNodes(); a++) {
       *pBdil++ = *pDbx*fthird; *pBdev++ = *pDx*2.0*fthird;
      *pBdil++ = *pDbx*fthird; *pBdev++ = *pDx*(-fthird);
      *pBdil++ = 0.0;          *pBdev++ = *pDy;
      *pBdil++ = *pDbx*fthird; *pBdev++ = *pDx*(-fthird);
      
      *pBdil++ = *pDby*fthird; *pBdev++ = *pDy*(-fthird);
      *pBdil++ = *pDby*fthird; *pBdev++ = *pDy*2.0*fthird;
      *pBdil++ = 0.0;          *pBdev++ = *pDx;
      *pBdil++ = *pDby*fthird; *pBdev++ = *pDy*(-fthird);
    
      pDx++;  pDbx++; 
      pDy++;  pDby++;
    }
  }
  else if (DNabar.MajorDim() == 3) {
    const double* pDx = DNa(0); const double* pDbx = DNabar(0); 
    const double* pDy = DNa(1); const double* pDby = DNabar(1);  
    const double* pDz = DNa(2); const double* pDbz = DNabar(2);

    for (int a = 0; a<NumElementNodes(); a++) {
      *pBdil++ = *pDbx*fthird; *pBdev++ = *pDx*2.0*fthird; 
      *pBdil++ = *pDbx*fthird; *pBdev++ = *pDx*(-fthird);
      *pBdil++ = *pDbx*fthird; *pBdev++ = *pDx*(-fthird);
      *pBdil++ = 0.0;          *pBdev++ = *pDz;
      *pBdil++ = 0.0;          *pBdev++ = 0.0;
      *pBdil++ = 0.0;          *pBdev++ = *pDy;
      
      *pBdil++ = *pDby*fthird; *pBdev++ = *pDy*(-fthird);
      *pBdil++ = *pDby*fthird; *pBdev++ = *pDy*2.0*fthird;
      *pBdil++ = *pDby*fthird; *pBdev++ = *pDy*(-fthird);
      *pBdil++ = 0.0;          *pBdev++ = 0.0;
      *pBdil++ = 0.0;          *pBdev++ = *pDz;
      *pBdil++ = 0.0;          *pBdev++ = *pDx;
      
      *pBdil++ = *pDbz*fthird; *pBdev++ = *pDz*(-fthird);
      *pBdil++ = *pDbz*fthird; *pBdev++ = *pDz*(-fthird);
      *pBdil++ = *pDbz*fthird; *pBdev++ = *pDy*(2.0*fthird);
      *pBdil++ = 0.0;          *pBdev++ = *pDx;
      *pBdil++ = 0.0;          *pBdev++ = *pDy;
      *pBdil++ = 0.0;          *pBdev++ = 0.0;
      
      pDx++;  pDbx++; 
      pDy++;  pDby++;
      pDz++;  pDbz++;
    }
  }
  else ExceptionT::OutOfRange("SmallStrainQ2P1::Set_Bbar");

  Bbar = B_dev;
  Bbar += Bbar_dil;
}

/* calculate the internal force contribution ("-k*d") */
void SmallStrainQ2P1::FormKd(double constK)
{
  fAVec = 0.0;  
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  
  fShapes->TopIP();
  while (fShapes->NextIP())  {
    Set_Bbar(fShapes->Derivatives_U(), fGradbar[CurrIP()],fBbar, fB_dev, fBbar_dil);
    const dSymMatrixT& stress = fCurrMaterial->s_ij();
 
   /*calculate pressure*/
    double pressure = fthird*stress.Trace();
    fMShapes.RowCopy(CurrIP(),fGamma);
    fAVec.AddScaled(pressure*(*Weight)*(*Det),fGamma);

    /*internal force*/
    fBbar.MultTx(stress, fNEEvec);    
    fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);		
  }	

  int elem = CurrElementNumber();  
  for (int ip = 0; ip<NumIP(); ip++)  {
    fMShapes.RowCopy(ip,fGamma);
    fPressure_List(elem,ip) = fH_inv.MultmBn(fGamma,fAVec);
  }
}

void SmallStrainQ2P1::FormStiffness(double constK)
{
  /* matrix format */
  dMatrixT::SymmetryFlagT format =
    (fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
    dMatrixT::kWhole :
    dMatrixT::kUpperOnly;
  
  /* integrate element stiffness */
  const double* Det    = fShapes->IPDets();
  const double* Weight = fShapes->IPWeights();
  
  fShapes->TopIP();
  while ( fShapes->NextIP() ) {
    double scale = constK*(*Det++)*(*Weight++);
    Set_Bbar(fShapes->Derivatives_U(), fGradbar[CurrIP()],fBbar, fB_dev, fBbar_dil);
    
    /* get D matrix */
    const dMatrixT& modulus = fCurrMaterial->c_ijkl();
    fD.SetToScaled(scale, modulus);
    
    /* multiply b(transpose) * db, taking account of symmetry, */
    /* and accumulate in elstif */
    fLHS.MultQTBQ(fBbar, fD, format, dMatrixT::kAccumulate);
  }
}

