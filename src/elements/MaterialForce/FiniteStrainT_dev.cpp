/* $Id: FiniteStrainT_dev.cpp,v 1.2 2003-02-12 00:36:50 thao Exp $ */
#include "FiniteStrainT.h"

#include "ScheduleT.h"
#include "Traction_CardT.h"
#include "ShapeFunctionT.h"
#include "FSSolidMatT.h"
#include "FSMatSupportT.h"
#include "ModelManagerT.h"

/* materials lists */
#include "SolidMatList1DT.h"
#include "SolidMatList2DT.h"
#include "SolidMatList3DT.h"

using namespace Tahoe;

/* constructor */
FiniteStrainT::FiniteStrainT(const ElementSupportT& support, const FieldT& field):
	SolidElementT(support, field),
	fNeedsOffset(-1),
	fCurrShapes(NULL),
	fFSMatSupport(NULL)
{
	/* disable any strain-displacement options */
	if (fStrainDispOpt != kStandardB)
	{
		cout << "\n FiniteStrainT::FiniteStrainT: no strain-displacement options\n" << endl;
		fStrainDispOpt = kStandardB;
	}
}

/* destructor */
FiniteStrainT::~FiniteStrainT(void)
{
	delete fFSMatSupport;
}

/* called immediately after constructor */
void FiniteStrainT::Initialize(void)
{
	/* inherited */
	SolidElementT::Initialize();

	/* what's needed */
	bool need_F = false;
	bool need_F_last = false;
	for (int i = 0; i < fMaterialList->Length(); i++)
	{
		need_F = need_F || Needs_F(i);		
		need_F_last = need_F_last || Needs_F_last(i);
	}	

	/* allocate deformation gradient list */
	if (need_F)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_all.Dimension(nip*nsd*nsd);
		fF_List.Dimension(nip);
		for (int i = 0; i < nip; i++)
			fF_List[i].Set(nsd, nsd, fF_all.Pointer(i*nsd*nsd));
	}
	
	/* allocate "last" deformation gradient list */
	if (need_F_last)
	{
		int nip = NumIP();
		int nsd = NumSD();
		fF_last_all.Dimension(nip*nsd*nsd);
		fF_last_List.Dimension(nip);
		for (int i = 0; i < nip; i++)
			fF_last_List[i].Set(nsd, nsd, fF_last_all.Pointer(i*nsd*nsd));
	}
}

/* TEMPORARY */
void FiniteStrainT::InitialCondition(void)
{
	/* inherited */
	SolidElementT::InitialCondition();
	
	/* set the source for the iteration number */
	fFSMatSupport->SetIterationNumber(ElementSupport().IterationNumber(Group()));
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u) const
{
	if (fCurrShapes)
	{
		/* field gradient */
		fCurrShapes->GradU(u, grad_u);
	}
	else
	{
		cout << "\n FiniteStrainT::ComputeGradient: shape functions wrt current coords not defined" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* compute field gradients with respect to current coordinates */
void FiniteStrainT::ComputeGradient(const LocalArrayT& u, dMatrixT& grad_u, 
	int ip) const
{
	if (fCurrShapes)
	{
		/* field gradient */
		fCurrShapes->GradU(u, grad_u, ip);
	}
	else
	{
		cout << "\n FiniteStrainT::ComputeGradient: shape functions wrt current coords not defined" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct a new material support and return a pointer */
MaterialSupportT* FiniteStrainT::NewMaterialSupport(MaterialSupportT* p) const
{
	/* allocate */
	if (!p) p = new FSMatSupportT(NumSD(), NumDOF(), NumIP());

	/* inherited initializations */
	SolidElementT::NewMaterialSupport(p);
	
	/* set FiniteStrainT fields */
	FSMatSupportT* ps = dynamic_cast<FSMatSupportT*>(p);
	if (ps) {
		ps->SetDeformationGradient(&fF_List);
		ps->SetDeformationGradient_last(&fF_last_List);
	}

	return p;
}

/* construct materials manager and read data */
MaterialListT* FiniteStrainT::NewMaterialList(int size)
{
	/* material support */
	if (!fFSMatSupport) {
		fFSMatSupport = dynamic_cast<FSMatSupportT*>(NewMaterialSupport());
		if (!fFSMatSupport) throw ExceptionT::kGeneralFail;
	}

	if (NumSD() == 1)
		return new SolidMatList1DT(size, *fFSMatSupport);
	else if (NumSD() == 2)
		return new SolidMatList2DT(size, *fFSMatSupport);
	else if (NumSD() == 3)
		return new SolidMatList3DT(size, *fFSMatSupport);
	else
		return NULL;			
}

/* construct list of materials from the input stream */
void FiniteStrainT::ReadMaterialData(ifstreamT& in)
{
	/* inherited */
	SolidElementT::ReadMaterialData(in);

	/* offset to class needs flags */
	fNeedsOffset = fMaterialNeeds[0].Length();
	
	/* set material needs */
	for (int i = 0; i < fMaterialNeeds.Length(); i++)
	{
		/* needs array */
		ArrayT<bool>& needs = fMaterialNeeds[i];

		/* resize array */
		needs.Resize(needs.Length() + 2, true);

		/* casts are safe since class contructs materials list */
		ContinuumMaterialT* pcont_mat = (*fMaterialList)[i];
		FSSolidMatT* mat = (FSSolidMatT*) pcont_mat;

		/* collect needs */
		needs[fNeedsOffset + kF     ] = mat->Need_F();
		needs[fNeedsOffset + kF_last] = mat->Need_F_last();
		
		/* consistency */
		needs[kNeedDisp] = needs[kNeedDisp] || needs[fNeedsOffset + kF];
		needs[KNeedLastDisp] = needs[KNeedLastDisp] || needs[fNeedsOffset + kF_last];
	}
}

/* form shape functions and derivatives */
void FiniteStrainT::SetGlobalShape(void)
{
	/* inherited */
	SolidElementT::SetGlobalShape();

	/* what needs to get computed */
	int material_number = CurrentElement().MaterialNumber();
	bool needs_F = Needs_F(material_number);
	bool needs_F_last = Needs_F_last(material_number);
	
	/* loop over integration points */
	for (int i = 0; i < NumIP(); i++)
	{
		/* deformation gradient */
		if (needs_F)
		{
			dMatrixT& mat = fF_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}

		/* "last" deformation gradient */
		if (needs_F_last)
		{
			dMatrixT& mat = fF_last_List[i];

			/* displacement gradient */
			fShapes->GradU(fLocLastDisp, mat, i);

			/* add identity */
			mat.PlusIdentity();
		}
	}
}

/* write all current element information to the stream */
void FiniteStrainT::CurrElementInfo(ostream& out) const
{
	/* inherited */
	SolidElementT::CurrElementInfo(out);
	
	/* write deformation gradients */
	out << "\n i.p. deformation gradients:\n";
	for (int i = 0; i < fF_List.Length(); i++)
		out << " ip: " << i+1 << '\n'
		    << fF_List[i] << '\n';
	out << '\n';
}


/*****************************************************************************/
/*public                                                                     */
/*****************************************************************************/
bool FiniteStrainT::MatForceDriver(void)
{
  /*obtain dimensions*/
  ModelManagerT& model = ElementSupport().Model();
  int nnd = model.NumNodes();
  int nen = NumElementNodes();
  int nsd = NumSD();
  int nmf = nnd*nsd;

  fMatForce.Dimension(nmf);
  fMatForceDissip.Dimension(nmf);
  fMatForce = 0;
  fMatForceDissip = 0;

  dArrayT MatForce(nsd*nen);
  int dissip = 0;
  
  Top();
  while (ContinuumElementT::NextElement())
  {
    ContinuumMaterialT* pmat = (*fMaterialList)[CurrentElement().MaterialNumber()];
    FSSolidMatT* CurrMaterial = (FSSolidMatT*)pmat; //Q:cast OK?
    /*Set Global Shape Functions for current element*/
    SetGlobalShape();

    MatForceVolMech(CurrMaterial, MatForce);
    AssembleMatForce(MatForce, fMatForce, CurrentElement().NodesX());
    if (CurrMaterial->HasDissipVar()) 
    {
      dissip++;
      MatForceDissip(CurrMaterial, MatForce, CurrentElement().DoubleData());
      AssembleMatForce(MatForce, fMatForceDissip, CurrentElement().NodesX());
    }
  }
  if (dissip >0)
  {
    fMatForce += fMatForceDissip;
    return(true);
  }
  else return(false);
}

void FiniteStrainT::AssembleMatForce(const dArrayT& elem_val, 
				     dArrayT& global_val, 
				     const nArrayT<int>& ndnos)
{
  int num_entries = global_val.Length();
  int nen = ndnos.Length();
  int nsd = elem_val.Length()/nen;
  int first_node = 1; //Q:assume first node number is 1
  for(int i = 0; i < nen; i++)
  {
    int node = ndnos[i];
    double* p = global_val.Pointer()+(node-first_node);
    for (int j = 0; j<nsd; j++)
      (*p++) += elem_val[i*nsd+j];
  }
}


void FiniteStrainT::MatForceVolMech(FSSolidMatT* CurrMaterial, 
				    dArrayT& MatForce)
{
  /*obtain dimensions*/
  int nsd = NumSD();
  int nen = NumElementNodes();

  /*initialize workspaces*/
  //element vectors: F=[F1x|F1y||F2x|F2y||...]
  dArrayT bodyforce(nsd*nen);
  dMatrixT Eshelby(nsd);
  dSymMatrixT C(nsd);
  dArrayT ip_body(nsd);
 
  MatForce = 0;
  
  /*get density*/
  double density = CurrMaterial->Density();
  
  /*obtain local nodal data*/
  SetLocalArrays();
  //Q:Is fLocAcc initialized to zero?
  if (fLocAcc.IsRegistered())
    SetLocalU(fLocAcc);

  /*copy acceleration and body force data into body force vector*/
  double* pbody = bodyforce.Pointer();    
  if (fBodySchedule)
  {
    double loadfactor = fBodySchedule->Value();
    for (int i = 0; i<nen; i++)
    {
      double* pLocAcc = fLocAcc.Pointer();
      pLocAcc += i;
      for (int j = 0; j<nsd; j++)
      {
	*pbody++ = (*pLocAcc - fBody[j]*loadfactor)*density;
	pLocAcc += nen;
      }
    }
  }
  else 
  {
    for (int i = 0; i<nen; i++)
    {
      double* pLocAcc = fLocAcc.Pointer();
      pLocAcc += i;
      for (int j = 0; j<nsd; j++)
      {
	*pbody++ = (*pLocAcc)*density;
	pLocAcc += nen;
      }
    }
  }
    
  /*intialize shape function data*/
  const double* jac = fShapes->IPDets();
  const double* weight = fShapes->IPWeights();
  fShapes->TopIP();
  while(fShapes->NextIP())
  {
    /*get shape function and derivatives at integration point*/
    const double* pQa = fShapes->IPShapeX();
    const double* pQaU = fShapes->IPShapeU();
    const dArray2DT& DQa = fShapes->Derivatives_X();
    
    /*gather material ip data*/
    double energy = CurrMaterial->StrainEnergyDensity();
    const dMatrixT& F = CurrMaterial->F_mechanical();
    Eshelby = 0;
    const dSymMatrixT& S = CurrMaterial->S_IJ();
    double* pbody = bodyforce.Pointer();
    double* pforce = MatForce.Pointer(); 
    if (nsd == 2)
    {
      /*interpolate to ip*/
      ip_body = 0;
      for (int i= 0; i<nen; i++)
      {
	ip_body[0] += (*pQaU) * (*pbody++);
	ip_body[1] += (*pQaU++) * (*pbody++);
      }	 
      /*form right stretch tensor*/
      C[0] = F[0]*F[0] + F[1]*F[1];
      C[1] = F[2]*F[2] + F[3]*F[3];
      C[2] = F[0]*F[2] + F[1]*F[3];
      /*form negative of Eshelby stress SIG_IJ = C_IK S_KJ - Psi Delta_IJ*/
      Eshelby[0] = C[0]*S[0]+C[2]*S[2]-energy; //(0,0) 
      Eshelby[1] = C[2]*S[0]+C[1]*S[2]; //(1,0)
      Eshelby[2] = C[0]*S[2]+C[2]*S[1]; //(0,1)
      Eshelby[3] = C[2]*S[2]+C[1]*S[1]-energy; //(1,1)
      
      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      
      for (int j = 0; j<nen; j++)
      {
	/*add Eshelby volume integral contribution*/
	*(pforce++) += ( Eshelby[0]*(*pDQaX)+Eshelby[2]*(*pDQaY)
		     + (F[0]*ip_body[0]+F[1]*ip_body[1])*(*pQa) )
	             *(*jac)*(*weight);
	*(pforce++) += ( Eshelby[1]*(*pDQaX++)+Eshelby[3]*(*pDQaY++)
		     + (F[2]*ip_body[0]+F[3]*ip_body[1])*(*pQa++) )
	             *(*jac)*(*weight);
      }
    }
    else if (nsd ==3)
    {
      /*interpolate ip values*/
      ip_body = 0;
      for (int i= 0; i<nen; i++)
      {
	ip_body[0] += (*pQaU) * (*pbody++);
	ip_body[1] += (*pQaU) * (*pbody++);
	ip_body[2] += (*pQaU++) * (*pbody++);
      }
	  
      /*form Eshelby stress SIG_IJ = Psi I_IJ - C_IK S_KJ*/
      C[0] = F[0]*F[0]+F[1]*F[1]+F[2]*F[2];
      C[1] = F[3]*F[3]+F[4]*F[4]+F[5]*F[5];
      C[2] = F[6]*F[6]+F[7]*F[7]+F[8]*F[8];
      C[3] = F[3]*F[6]+F[4]*F[7]+F[5]*F[8];
      C[4] = F[0]*F[6]+F[1]*F[7]+F[2]*F[8];
      C[5] = F[0]*F[3]+F[1]*F[4]+F[2]*F[5];
      
      Eshelby[0] = C[0]*S[0]+C[5]*S[5]+C[4]*S[4]-energy; //(0,0)
      Eshelby[4] = C[5]*S[5]+C[1]*S[1]+C[3]*S[3]-energy; //(1,1)
      Eshelby[8] = C[4]*S[4]+C[3]*S[3]*C[2]*S[2]-energy; //(2,2)
      Eshelby[7] = C[5]*S[4]+C[1]*S[3]+C[3]*S[2]; //(1,2)
      Eshelby[5] = C[4]*S[5]+C[3]*S[1]+C[2]*S[3]; //(2,1)
      Eshelby[6] = C[0]*S[4]+C[5]*S[3]+C[4]*S[2]; //(0,2)
      Eshelby[2] = C[4]*S[0]+C[3]*S[5]+C[2]*S[4]; //(2,0)
      Eshelby[3] = C[0]*S[5]+C[5]*S[1]*C[4]*S[3]; //(0,1)
      Eshelby[1] = C[5]*S[0]+C[1]*S[5]+C[3]*S[4]; //(1,0)
      
      double* pDQaX = DQa(0); 
      double* pDQaY = DQa(1);
      double* pDQaZ = DQa(2);
      
      for (int j = 0; j<nen; j++)
      {
	/*add Eshelby volume integral contribution*/
	*(pforce++) += (Eshelby[0]*(*pDQaX)+Eshelby[3]*(*pDQaY)
	             +  Eshelby[6]*(*pDQaZ)
                     + (F[0]*ip_body[0]+F[1]*ip_body[1]
                     +  F[2]*ip_body[2])*(*pQa) )*(*jac)*(*weight);
	*(pforce++) += (Eshelby[1]*(*pDQaX)+Eshelby[4]*(*pDQaY)
                     +  Eshelby[7]*(*pDQaZ)
                     + (F[3]*ip_body[0]+F[4]*ip_body[1]
		     +	F[5]*ip_body[2])*(*pQa) )*(*jac)*(*weight);
	*(pforce++) += (Eshelby[2]*(*pDQaX++)+Eshelby[5]*(*pDQaY++)
		     +  Eshelby[8]*(*pDQaZ++)
		     + (F[6]*ip_body[0]+F[7]*ip_body[1]
                     +  F[8]*ip_body[2])*(*pQa++))*(*jac)*(*weight);
      }
    }
    weight++;
    jac++;
  }
}

void FiniteStrainT::MatForceDissip(FSSolidMatT* CurrMaterial,dArrayT& MatForce,
				   dArrayT& statev)
{
  /*obtain dimensions*/
  int nsd = NumSD();
  int numstress = nsd*(nsd+1)*0.5;
  int nen = NumElementNodes();
  
  /*initialize workspaces*/
  dSymMatrixT ExtrapMatrix(nen);
  nArrayT<dSymMatrixT> Nodal_iInStretch(nen);
  nArrayT<dSymMatrixT> IP_iInStretch(nen);
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& ip_val = IP_iInStretch[i];
    dSymMatrixT& nodal_val = Nodal_iInStretch[i];
    ip_val.Allocate(nsd);
    nodal_val.Allocate(nsd);
  }
  dArrayT Grad_iInStretch(nsd*numstress);

  MatForce = 0;
    
  /*get shape function data*/
  int nip = fShapes->NumIP();
  const double* jac = fShapes->IPDets();;
  const double* weight = fShapes->IPWeights();

  /*set inelastic stretch to current values*/
  int nstatev = statev.Length()/nip;
  double* pstatev = statev.Pointer();
  double* pInStretch = pstatev;
  pstatev += numstress;
  /*skip over previous value of inelastic stretch
    and set inelastic stress measure*/
  pstatev += numstress;
  double* pInStress = pstatev;

  /*extrapolate iInStretch to nodes*/
  ExtrapMatrix = 0;
  /*get shape function and derivatives at integration point*/
  fShapes->TopIP();
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& ip_val = IP_iInStretch[i];
    ip_val = 0;
  }
  while(fShapes->NextIP())
  {
    dSymMatrixT InStretch(numstress, pInStretch);  //Q:Is this sufficient?
    dSymMatrixT iInStretch = InStretch;
    iInStretch.Inverse();

    const double* pQbU = fShapes->IPShapeU();
    for (int i=0; i<nen; i++)
    {
      const double* pQaU = fShapes->IPShapeU();
      for (int j = i; j<nen; i++)
      {
	//Note: The extrapolation matrix is symmetric
	ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
      }
      
      dSymMatrixT& ip_val = IP_iInStretch[i];
      ip_val.AddCombination((*pQbU++)*(*jac)*(*weight), iInStretch);
    }
    pInStretch += nstatev;
    weight++;
    jac++;
  }
  ExtrapMatrix.Inverse();
  for (int i = 0; i<nen; i++)
  {
    dSymMatrixT& nodal_val = Nodal_iInStretch[i];
    nodal_val = 0;
    for (int j = 0; j<nen; j++)
      nodal_val.AddCombination(ExtrapMatrix[i+j*nen],IP_iInStretch[j]);
  }
  
  double* pMatForce = MatForce.Pointer();
  jac = fShapes->IPDets();
  weight = fShapes->IPWeights();
  fShapes->TopIP();
  while(fShapes->NextIP())
  {
    const double* pQa = fShapes->IPShapeX();
    const dArray2DT& DQa = fShapes->Derivatives_X();
    if (nsd ==2)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);
      Grad_iInStretch = 0;

      /*Interpolate grad of iverse inelastic stretch to ip*/
      double* pGradX = Grad_iInStretch.Pointer();
      double* pGradY = pGradX+numstress;
      for (int i = 0; i<nen; i++)
      {
	dSymMatrixT& nodal_val = Nodal_iInStretch[i];

        pGradX[0] += (*pDQaX)*nodal_val[0];
	pGradX[1] += (*pDQaX)*nodal_val[1];
	pGradX[2] += (*pDQaX++)*nodal_val[2];
	
	pGradY[0] += (*pDQaY)*nodal_val[0];
	pGradY[1] += (*pDQaY)*nodal_val[1];
	pGradY[2] += (*pDQaY++)*nodal_val[2];
      }
      /*get inelastic stress*/
      dSymMatrixT InStress(numstress, pInStress);
      
      /*integrate material force*/
      for (int i = 0; i<nen; i++)
      {
	(*pMatForce++) +=0.5*(InStress[0]*pGradX[0]+InStress[1]*pGradX[1]
			+2.0* InStress[2]*pGradX[2])*(*pQa)*(*jac)*(*weight);
	(*pMatForce++) +=0.5*(InStress[0]*pGradY[0]+InStress[1]*pGradY[1]
		        +2.0* InStress[2]*pGradY[2])*(*pQa++)*(*jac)*(*weight);
      }
    }
    else if (nsd ==3)
    {
      const double* pDQaX = DQa(0);
      const double* pDQaY = DQa(1);
      const double* pDQaZ = DQa(2);
      Grad_iInStretch = 0;
      /*Interpolate grad of iverse inelastic stretch to ip*/
      double* pGradX = Grad_iInStretch.Pointer();
      double* pGradY = pGradX+numstress;
      double* pGradZ = pGradY+numstress;
      for (int i = 0; i<nen; i++)
      {	
	dSymMatrixT& nodal_val = Nodal_iInStretch[i];

	pGradX[0] += (*pDQaX)*nodal_val[0];
	pGradX[1] += (*pDQaX)*nodal_val[1];
	pGradX[2] += (*pDQaX)*nodal_val[2];
	pGradX[3] += (*pDQaX)*nodal_val[3];
	pGradX[4] += (*pDQaX)*nodal_val[4];
	pGradX[5] += (*pDQaX++)*nodal_val[5];
	
	pGradY[0] += (*pDQaY)*nodal_val[0];
	pGradY[1] += (*pDQaY)*nodal_val[1];
	pGradY[2] += (*pDQaY)*nodal_val[2];
	pGradY[3] += (*pDQaY)*nodal_val[3];
	pGradY[4] += (*pDQaY)*nodal_val[4];
	pGradY[5] += (*pDQaY++)*nodal_val[5];
      }
      /*get inelastic stress*/
      dSymMatrixT InStress(numstress, pInStress);
      
      /*integrate material force*/
      for (int i = 0; i<nen; i++)
      {
	(*pMatForce++) += 0.5*(InStress[0]*pGradX[0]+InStress[1]*pGradX[1]
 		           +   InStress[2]*pGradX[2]+2.0*InStress[3]*pGradX[3]
                          +2.0*InStress[4]*pGradX[4]+2.0*InStress[5]*pGradX[5])
                           *(*pQa)*(*jac)*(*weight);
	(*pMatForce++) += 0.5*(InStress[0]*pGradY[0]+InStress[1]*pGradY[1]
	                   +   InStress[2]*pGradY[2]+2.0*InStress[3]*pGradY[3]
                          +2.0*InStress[4]*pGradY[4]+2.0*InStress[5]*pGradY[5])
                           *(*pQa)*(*jac)*(*weight);
	(*pMatForce++) += 0.5*(InStress[0]*pGradZ[0]+InStress[1]*pGradZ[1]
	                   +   InStress[2]*pGradZ[2]+2.0*InStress[3]*pGradZ[3]
                          +2.0*InStress[4]*pGradZ[4]+2.0*InStress[5]*pGradZ[5])
                           *(*pQa++)*(*jac)*(*weight);
      }
    }
    jac++;
    weight++;
  }
}

void FiniteStrainT::MatForceSurfMech(void)
{
  if(fTractionList.Length() > 0)
  {
    /*obtain dimensions*/
    int nen = NumElementNodes();
    int nsd = NumSD();

    /*initialize workspace*/
    dSymMatrixT ExtrapMatrix(nen);
    nArrayT<dMatrixT> IP_F(nen);
    for (int i = 0; i<nen; i++)
    {
      dMatrixT& ip_val = IP_F[i];
      ip_val.Allocate(nsd);
    }

    dArrayT IP_energy(nen);      

    dMatrixT jacobian(nsd,nsd-1);
    dMatrixT ip_F(nsd);
    dArrayT ip_tract(nsd);
    dArrayT ip_mattract(nsd);
    dArrayT global_tract(nsd);
    double ip_energy;
    dMatrixT Q(nsd);
 
    /*update eqnos*/
    //Q:need to set tractions: i.e. SetTractionBC()?
  
    for (int k = 0; k<fTractionList.Length(); k++)
    {
      /*retrieve ith traction card*/
      const Traction_CardT& BC_card = fTractionList[k];

      /*dimension of ith traction card*/
      const iArrayT& surf_nodes = BC_card.Nodes(); //Q:global node numbers
      int nsn = surf_nodes.Length();
      /*LocalArrayT: nenxnsd F=[F1X|F2X|F3X|...||F1Y|F2Y|F3Y|....]*/
      /*retrieve coordinates of nodes in local ordering*/
      LocalArrayT surf_coords(LocalArrayT::kInitCoords, nsn, nsd);
      ElementSupport().RegisterCoordinates(surf_coords);
      surf_coords.SetLocal(surf_nodes);    

      /*allocate work space*/
      nArrayT<dMatrixT> Nodal_F(nsn);
      for (int i = 0; i<nsn; i++)
      {
	dMatrixT& nodal_val = Nodal_F[i];
	nodal_val.Allocate(nsd);
      }
      dArrayT nodal_energy(nsn);  

      /*initialize mat force vector*/
      dArrayT MatForce(nsn*nsd);
      /*get local nodal tractions scaled by LTf*/
      LocalArrayT tract(LocalArrayT::kUnspecified, nsn, nsd);
      BC_card.CurrentValue(tract);

      /*BC destination*/
      int elem, facet;  //element # and facet # of traction bc
      BC_card.Destination(elem,facet);

      /*Extrapolate ip energies and def grad to nodes*/
      
      const ElementCardT& elem_card = fElementCards[elem];
      ContinuumMaterialT* pmat = (*fMaterialList)[elem_card.MaterialNumber()];
      FSSolidMatT* CurrMaterial = (FSSolidMatT*) pmat; //Q: is this right?

      /*get element coords*/
      const iArrayT& elem_nodes=elem_card.NodesX();
      LocalArrayT elem_coords(LocalArrayT::kInitCoords,nen,nsd);
      ElementSupport().RegisterCoordinates(elem_coords);
      elem_coords.SetLocal(elem_nodes);    

      /*Set element shapefunctions*/    //Q:How do we do this? 
      fShapes->SetDerivatives(); 
      const double* jac = fShapes->IPDets();
      const double* weight = fShapes->IPWeights();

      /*Form Extrapolation Matrix*/
      ExtrapMatrix = 0;
      IP_energy = 0;
      for (int i = 0; i<nen; i++)
      {
	dMatrixT& ip_val = IP_F[i];
	ip_val = 0;
      }
      fShapes->TopIP();
      while (fShapes->NextIP())
      {
	const dMatrixT& F = CurrMaterial->F_mechanical();
	double energy = CurrMaterial->StrainEnergyDensity();

	const double* pQbU = fShapes->IPShapeU();
	for (int i=0; i<nen; i++)
	{
	  const double* pQaU = fShapes->IPShapeU();
	  for (int j = i; j<nen; i++)
	  {
	    //Note: The extrapolation matrix is symmetric
	    ExtrapMatrix(i,j) += (*pQaU++)*(*pQbU)*(*jac)*(*weight);  
	  }

	  IP_energy[i] += (*pQbU)*energy*(*jac)*(*weight);
	  dMatrixT& ip_val = IP_F[i];
	  ip_val.AddCombination((*pQbU++)*(*jac)*(*weight), F);
	}
	weight++;
	jac++;
      }
      ExtrapMatrix.Inverse();
      //Q:Assumes that nodes are in local ordering
      for (int i = 0; i<nsn; i++)
      {
	int cnt = 0;
	while(elem_nodes[cnt]!=surf_nodes[i]) cnt++; 
	dMatrixT& nodal_val = Nodal_F[i];
	nodal_val = 0;
	for (int j = 0; j<nen; j++)
	{
	  nodal_val.AddCombination(ExtrapMatrix[cnt+j*nen],IP_F[j]);
	  nodal_energy[i] += ExtrapMatrix[cnt+j*nen]*IP_energy[j];
	}
      }

      /*get surface shape function*/
      const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
      int nip = surf_shape.NumIP();       
      double* penerg = nodal_energy.Pointer();
       
      if (nsd ==2)
      {
	MatForce = 0.0;
	int thickness = 1.0; //Q:Hardwired thickness;
	const double* ip_w = surf_shape.Weight();
	for (int j = 0; j < nip; j++)
	{
	  surf_shape.DomainJacobian(surf_coords,j,jacobian);
	  double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	  
	  /*interpolate to surface ip*/
	  double* ptract_X = tract(0);
	  double* ptract_Y = tract(1);
	  ip_F = 0;
	  ip_tract = 0;
	  ip_energy = 0;
	  
	  const double* pQaU = surf_shape.Shape(j);
	  for (int i= 0; i<nsn; i++)
	  {
	    dMatrixT& nodal_val = Nodal_F[i];
	    ip_F.AddCombination(*pQaU,nodal_val);
	    ip_energy += (*pQaU)*(*penerg++);
	    ip_tract[0] += (*pQaU)*(*ptract_X++);
	    ip_tract[1] += (*pQaU++)*(*ptract_Y++);
	  }
	  /*surface normal*/
	  double* n = Q(1);
	  if (BC_card.CoordSystem() == Traction_CardT::kLocal)
	  {
	    /*rotate traction from local to global coords*/
	    global_tract[0] = Q[0]*ip_tract[0]+Q[2]*ip_tract[1];
	    global_tract[1] = Q[1]*ip_tract[0]+Q[3]*ip_tract[1];
	  }
	  else if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    global_tract[0] = ip_tract[0];
	    global_tract[1] = ip_tract[1];
	  }
	  ip_mattract[0] = (ip_F[0]*global_tract[0]
			    +ip_F[1]*global_tract[1])*thickness-ip_energy*n[0];
	  ip_mattract[1] = (ip_F[2]*global_tract[0]
			    +ip_F[3]*global_tract[1])*thickness-ip_energy*n[1];
	  
	  /*integrate material force*/
	  const double* pQa = surf_shape.Shape(j);
	  double * pMatForce = MatForce.Pointer();
	  for (int i = 0; i<nsn; i++)
	  {
	    (*pMatForce++) +=(*pQa)*ip_mattract[0]*detj*(*ip_w);
	    (*pMatForce++) +=(*pQa++)*ip_mattract[1]*detj*(*ip_w);
	  }
	  ip_w++;
	}
      }
      else if (nsd == 3)
      {
	const double* ip_w = surf_shape.Weight();
	for (int j = 0; j < nip; j++)
	{
	  surf_shape.DomainJacobian(surf_coords,j,jacobian);
	  double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	  
	  /*interpolate to surface ip*/
	  double* ptract_X = tract(0);
	  double* ptract_Y = tract(1);
	  double* ptract_Z = tract(2);
	  ip_F = 0;
	  ip_tract = 0;
	  ip_energy = 0;
	  const double* pQaU = surf_shape.Shape(j);
	  for (int i= 0; i<nsn; i++)
	  {
	    dMatrixT& nodal_val = Nodal_F[i];
	    ip_F.AddCombination(*pQaU,nodal_val);
	    ip_energy += (*pQaU)*(*penerg++);
	    
	    ip_tract[0] += (*pQaU)*(*ptract_X++);
	    ip_tract[1] += (*pQaU)*(*ptract_Y++);
	    ip_tract[2] += (*pQaU++)*(*ptract_Z++);
	  }
	  /*surface normal*/
	  double* n = Q(2);
	  if (BC_card.CoordSystem() == Traction_CardT::kLocal)
	  {
	    /*rotate traction from local to global coords*/
	    global_tract[0]=Q[0]*ip_tract[0]+Q[3]*ip_tract[1]+Q[6]*ip_tract[2];
	    global_tract[1]=Q[1]*ip_tract[0]+Q[4]*ip_tract[1]+Q[7]*ip_tract[2];
	    global_tract[2]=Q[2]*ip_tract[0]+Q[5]*ip_tract[1]+Q[8]*ip_tract[2];
	  }
	  else if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
	  {
	    global_tract[0] = ip_tract[0];
	    global_tract[1] = ip_tract[1];
	    global_tract[2] = ip_tract[2];
	  }
	  
	  ip_mattract[0] = (ip_F[0]*ip_tract[0]+ip_F[1]*ip_tract[1]
			  + ip_F[2]*ip_tract[2]) - ip_energy*n[0];
	  ip_mattract[1] = (ip_F[3]*ip_tract[0]+ip_F[4]*ip_tract[1]
			  + ip_F[5]*ip_tract[2]) - ip_energy*n[1];
	  ip_mattract[1] = (ip_F[6]*ip_tract[0]+ip_F[7]*ip_tract[1]
			  + ip_F[8]*ip_tract[2]) - ip_energy*n[2];

	  /*integrate material force*/
	  const double* pQa = surf_shape.Shape(j);
	  double * pMatForce = MatForce.Pointer();
	  for (int i = 0; i<nsn; i++)
	  {
	    (*pMatForce++) +=(*pQa)*ip_mattract[0]*detj*(*ip_w);
	    (*pMatForce++) +=(*pQa)*ip_mattract[1]*detj*(*ip_w);
	    (*pMatForce++) +=(*pQa++)*ip_mattract[2]*detj*(*ip_w);
	  }
	  ip_w++;
	}
      }
      //AssembleMatForce(Group(), MatForce, BC_card.Equations());
    }
  }
}
