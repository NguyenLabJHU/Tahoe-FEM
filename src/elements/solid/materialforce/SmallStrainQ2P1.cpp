/* $Id: SmallStrainQ2P1.cpp,v 1.1 2003-06-28 00:40:29 thao Exp $ */
#include "SmallStrainQ2P1.h"
#include "ShapeFunctionT.h"
#include "SSSolidMatT.h"

/* materials lists */
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"
#include "SSMatSupportT.h"

using namespace Tahoe;

/* constructor */
SmallStrainQ2P1::SmallStrainQ2P1(const ElementSupportT& support, const FieldT& field):
	SmallStrainT(support, field),
	fpdof(NumSD()+1),
	fStress(NumSD()),
	fModulus(dSymMatrixT::NumValues(NumSD())),
	fPdev(dSymMatrixT::NumValues(NumSD())),
	fOne(dSymMatrixT::NumValues(NumSD())),XS
	fthird(1.0/3.0)
{
        fMEArray.Dimension(fpdof);
	fMEArray = 0.0;
}

/* called immediately after constructor */
void SmallStrainQ2P1::Initialize(void)
{
	/* inherited */
	SmallStrainT::Initialize();
	
	int nsd = NumSD();
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
	     cout << "\n SimoQ1P0::Initialize: expecting quad geometry: "
		  << GeometryCode() << endl;
	     throw ExceptionT::kBadInputValue;
        }	

	/*initializes utility tensors*/
	if (nsd != 2){
	  cout << "\n SmallStrainQ2P1::Initialize: for now expects 2D plane strain geometry.";
	  throw ExceptionT::kBGeneralFail;
	}
	else {
	  fOne[0] = fOne[1] = 1.0; fOne[2] = 0.0;

	  fPdev = 0.0
	  fPdev(0,0) = fPdev(1,1) = 2.0*fthird;
	  fPdev(0,1) = fPdev(1,0) = fthird;
	  fPdev(2,2) = 1.0;
	}

	fTheta_List.Dimension(nel,nip);
	fTheta_List_last.Dimension(nel,nip);
	fPressure_List.Dimension(nel,nip);

	fTheta_List = 0.0;
	fTheta_List_last = 0.0;
	fPressure_List = 0.0;

	fMShapes.Dimension(nel);
	fH_inv.Dimension(nel);
	for (int i = 0; i < nel; i++)
	{
	     fMShapes[i].Dimension(nip, fpdof);
	     fH_inv[i].Dimension(fpdof);
	}
	fH_i.Dimension(fpdof);
	fbabar.Dimensions(nsd,nen);
	fGradNa.Dimension(nsd);
	fMixMat.Dimension(fpdof,nsd);
	fMixMat = 0.0;
}

void SmallStrainQ2P1::CloseStep(void)
{
        /*inherited*/
        ContinuumElementT::CloseStep();
	fTheta_List_last = fTheta_List;
}

void SmallStrainQ2P1::ResetStep(void)
{
        /*inherited*/
        ContinuumElementT::ResetStep();
	fTheta_List = fTheta_List_last;
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

/* calculate the internal force contribution ("-k*d") */
void SmallStrainQ2P1::FormKd(double constK)
{
        int elem = CurrElementNumber();

	const double* Det    = fShapes->IPDets();
	const double* Weight = fShapes->IPWeights();
	
	/********DEBUG*******/
	bool print = false; 
	int pos = fElementCards.Position(); 
	if (pos == 1&&0)  
	  print = true; 
	/*******************/
	
	/* mixed pressure*/
	fShapes->TopIP();
	fMEArray = 0.0;
	
	while (fShapes->NextIP())
	{
		/* Deviatoric Cauchy stress */
		double p = fCurrMaterial->mPressure();

		fMShapes[elem].RowCopy(CurrIP(),fGamma);
		fMEArray.AddScaled(p*(*Weight++)*(*Det++),fGamma);
	}	
		
	Det = fShapes->IPDets();
	Weight = fShapes->IPWeights();	
	double* ppres = fPressure_List(elem);
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
	        Set_B(fShapes->Derivatives_U(), fB);
		fGamma.Set(fpdof, fMShapes[elem](CurrIP()));

		fStress = fCurrMaterial->s_ij();
		double p = fCurrMaterial->mPressure();
		*ppres = fH_inv[elem].MultmBn(fGamma,fMEArray);;

		fStress.PlusIdentity((*ppres++)-p);
		fB.MultTx(fStress, fNEEvec);

		if (print) cout << "\nstress: "<<fStress;
		
		/* accumulate */
		fRHS.AddScaled(constK*(*Weight++)*(*Det++), fNEEvec);			}	
}

void SmallStrainQ2P1::FormStiffness(double constK)
{
        /* matrix format */
        dMatrixT::SymmetryFlagT format =
                (fLHS.Format() == ElementMatrixT::kNonSymmetric) ?
                dMatrixT::kWhole :
                dMatrixT::kUpperOnly;

        /********DEBUG*******/
        bool print = false; 
        int pos = fElementCards.Position(); 
        if (pos == 1&&0)  
          print = true; 
        /*******************/
	
	int nen = NumElementNodes();
	int elem = CurrElementNumber();

	for (int i = 0; i < nen; i++)
	{
	     /*set babar functions*/		    
	     const double* Det    = fShapes->IPDets();
	     const double* Weight = fShapes->IPWeights();

	     fShapes->TopIP();
	     while ( fShapes->NextIP() )
	     {
	           const dArray2DT& DQa = fShapes->Derivatives_U();
		   fMShapes[elem].RowCopy(CurrIP(),fGamma);
		   DQa.ColumCopy(i, fGradNa);
		   fMixMat.Outer(fGamma,fGradNa,(*Det++)*(*Weight++),nMatrixT::kAccumulate); 
	     }
	     fMixMat.MultAB(fH_inv[elem], fMixMat, nMatrixT::kWhole);
	}

        /* integrate element stiffness */
        Det    = fShapes->IPDets();
        Weight = fShapes->IPWeights();

        fShapes->TopIP();
        while ( fShapes->NextIP() )
        {

                double scale = constK*(*Det++)*(*Weight++);
        
		Set_B(fShapes->Derivatives_U(), fB);

                /* get D matrix */
		fModulus = fCurrMaterial->c_ijkl();
		
                fD.SetToScaled(scale, fCurrMaterial->c_ijkl());

                if (print) cout << "\nmodulus: "<<fCurrMaterial->c_ijkl();
                                                        
                /* multiply b(transpose) * db, taking account of symmetry, */
                /* and accumulate in elstif */
                fLHS.MultQTBQ(fB, fD, format, dMatrixT::kAccumulate);   
        }
 
} 

/* compute the measures of strain/deformation over the element */
void SmallStrainQ2P1::SetGlobalShape(void)
{
	/* inherited */
	SolidElementT::SetGlobalShape();

	/* material information */
	int material_number = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[material_number];

	/* element information */
	int nen = NumElementNodes();
	int nsd = NumSD();
	int nip = fShapes->NumIP();
	int elem = CurrElementNumber();

	const double* jac = fShapes->IPDets();
	const double* W = fShapes->IPWeights();  
	fH_inv[elem] = 0.0;
	fShapes->TopIP();
	while (fShapes->NextIP())
	{
	       int ip = CurrIP();

	       /* pressure shape functions */
	       fGamma.Set(fpdof, fMShapes[elem](ip));
	       fip_coords.Set(nsd,fGamma.Pointer(1));
	       fGamma[0] = 1.0;
	       IP_Coords(fip_coords);
	       
	       /* mass matrix formed from pressure shape functions */
	       fH_i.Outer(fGamma);
	       fH_inv[elem].AddScaled((*W)*(*jac),fH_i);
	       
	       /* deformation gradient */
	       if (needs[fNeedsOffset + kstrain])
	       {
		      /* displacement gradient */
		      fShapes->GradU(fLocDisp, fGradU, ip);

		      /* symmetric part */
		      fStrain_List[ip].Symmetrize(fGradU);
	       }

	       /* "last" deformation gradient */
	       if (needs[fNeedsOffset + kstrain_last])
	       {
		      /* displacement gradient */
		      fShapes->GradU(fLocLastDisp, fGradU, ip);

		      /* symmetric part */
		      fStrain_last_List[ip].Symmetrize(fGradU);
	       }      

	       /* evaluate dilation field*/
	       double dil = fStrain_List[ip].Trace();
	       for (int i = 0; i<fpdof; i++)
		      fMEArray[i] += fGamma[i]*dil*(*W++)*(*jac++);
	}
	fH_inv[elem].Inverse();

	/* set mixed strains */
	double* ptheta = fTheta_List(elem);
	double* ptheta_last = fTheta_List_last(elem);
	for(int i = 0; i < nip; i++)
	{	   
	       /* dilation field*/
	       fMShapes[elem].RowCopy(i,fGamma);
	       *ptheta = fH_inv[elem].MultmBn(fGamma,fMEArray);

	       /* deformation gradient */
	       if (needs[fNeedsOffset + kstrain])
	       {
		      /* correct for mixed field */
		      double dil = fStrain_List[i].Trace();
		      fStrain_List[i].PlusIdentity(fthird*( (*ptheta++)-dil));
	       }

	       /* "last" deformation gradient */
	       if (needs[fNeedsOffset + kstrain_last])
	       {
		      /*correct for mixed field*/
		      double dil = fStrain_last_List[i].Trace();
		      fStrain_last_List[i].PlusIdentity(fthird*( (*ptheta_last++)-dil));
	       }
	}
	      
}

