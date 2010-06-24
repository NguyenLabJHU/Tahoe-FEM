/* $Id: AnisoCornea.cpp,v 1.12 2010-06-24 13:49:16 thao Exp $ */
/* created: paklein (11/08/1997) */

#include "AnisoCornea.h"
#include <math.h>
#include "bessel.h"
#include "toolboxConstants.h"
#include "C1FunctionT.h"
#include "ParameterContainerT.h"
#include "ofstreamT.h"

#if defined(VIB_MATERIAL)

/* point generator */
#include "EvenSpacePtsT.h"

#include "FungType.h"
#include "PowerTrig.h"

const double Pi = acos(-1.0);
const double third = 1.0/3.0;
const int kNumOutputVar = 10;
static const char* Labels[kNumOutputVar] = 
	{"NT_X", "NT_Y", "NT_Z","IS_X", "IS_Y", "IS_Z","OT_X", "OT_Y", "OT_Z","J"};
static const int perm[3][3] = {0,1,2,1,2,0,2,0,1};

using namespace Tahoe;

/* constructors */
AnisoCornea::AnisoCornea(void):
    FSFiberMatSplitT(),
	ParameterInterfaceT("aniso_cornea"),
	fCircle(NULL),
	fPotential(NULL)
{
#ifndef VIB_MATERIAL
	ExceptionT::BadInputValue("AnisoCornea::AnisoCornea", 
		"VIB_MATERIAL must be enabled");
#endif

}

/* destructor */
AnisoCornea::~AnisoCornea(void) 
{ 
	delete fCircle; 
	delete fPotential;
}

/* strain energy density */
double AnisoCornea::StrainEnergyDensity(void)
{

	/*matrix contribution*/
	/* stretch */
	Compute_C(fC);
	double I3 = fC.Det();

	double I3rthird = pow(I3,-third);
	double I1bar = (fC[0]+fC[1]+fC[2]);
	I1bar *= I3rthird;
	
	/*coupled compressible Neo-Hookean*/
	/* mu/2 (I1 -3) +mu/(2 gamma) (I3^-gamma -1)*/
	double energy = 0.5*fMu*(I1bar-3.0) + 0.25*fKappa*(I3-1-log(I3));
	
	/*fiber contribution*/
	/* stretched bonds */
	ComputeLengths(fFiberStretch);

	/*equilibrium contribution*/
	/* update potential table */
	fPotential->MapFunction(fI4,fU);

	/* sum contributions */
	double* pU = fU.Pointer();
	double* pj =  fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

	for (int i = 0; i < fI4.Length(); i++)
		energy += (*pU++)*(*pj++);
	
	return energy;
}

int AnisoCornea::NumOutputVariables() const {
	return kNumOutputVar;
}

void AnisoCornea::OutputLabels(ArrayT<StringT>& labels) const
{
	//allocates space for labels
	labels.Dimension(kNumOutputVar);
	
	//copy labels
	for (int i = 0; i< kNumOutputVar; i++)
		labels[i] = Labels[i];
}

void AnisoCornea::ComputeOutput(dArrayT& output)
{
	/*calculates deformed fiber vectors*/
//	const dArray2DT& Fibers = FiberMatSupportT().Fiber_Vec();
	const dMatrixT& Q = GetRotation();
	const double* p_nt = Q(0);
	const double* p_is = Q(1);
	const double* p_ot = Q(2);
	
	const dMatrixT& F = F_mechanical();
	double* pb = output.Pointer();
	
	/*deformed NT fiber orientation*/
	F.Multx(p_nt, pb);
	pb += NumSD();
	
	/*deformed IS fiber orientation*/
	F.Multx(p_is, pb);
	pb += NumSD();

	/*deformed OT fiber orientation*/
	F.Multx(p_ot, pb);
	pb += NumSD();

	/*non-equilibrium strain energy density */
//	*pb = NonequilibriumStrainEnergyDensity();

	*pb = F.Det();
}

/* describe the parameters needed by the interface */
void AnisoCornea::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FSFiberMatSplitT::DefineParameters(list);
	
	/* integration points */
	ParameterT points(ParameterT::Integer, "n_points");
	points.AddLimit(1, LimitT::LowerInclusive);
	list.AddParameter(points);
}

/* information about subordinate parameter lists */
void AnisoCornea::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FSSolidMatT::DefineSubs(sub_list);

	/*material parameters for matrix*/
	sub_list.AddSub("matrix_material_params", ParameterListT::Once);

	/* choice of energy potential for fibrils */
	sub_list.AddSub("fibril_potential", ParameterListT::Once);

	/* choice of fibril distribution funcion */
	sub_list.AddSub("fibril_distribution", ParameterListT::Once);
}


/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* AnisoCornea::NewSub(const StringT& name) const
{
	/* inherited */
	ParameterInterfaceT* sub = FSSolidMatT::NewSub(name);

	if (sub) 
	{
		return sub;
	}
	else if (name == "matrix_material_params")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		/* bound */
		LimitT lower(0.0, LimitT::Lower);
		
		/* exponential functions*/
		ParameterContainerT matrix("Neo-Hookean");
		ParameterT mu(ParameterT::Double, "shear_modulus");
		mu.AddLimit(lower);
		ParameterT gamma(ParameterT::Double, "bulk_modulus");
		gamma.AddLimit(lower);
		matrix.AddParameter(mu);
		matrix.AddParameter(gamma);
		choice->AddSub(matrix);
		return choice;
	}
	else if (name == "fibril_potential")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);

		ParameterContainerT fung("fung_type");		
		{
			LimitT lower(0.0, LimitT::Lower);

			ParameterT alpha(ParameterT::Double, "alpha");
			ParameterT beta(ParameterT::Double, "beta");

			fung.AddParameter(alpha);
			fung.AddParameter(beta);
			alpha.AddLimit(lower);
			beta.AddLimit(lower);

			/* set the description */
			choice->AddSub(fung);
		}

		return(choice);
	}
	else if (name == "fibril_distribution")
	{
		ParameterContainerT* choice = new ParameterContainerT(name);
		choice->SetListOrder(ParameterListT::Choice);
		LimitT lower(0.0, LimitT::LowerInclusive);
		LimitT upper(1.0, LimitT::UpperInclusive);

		/* sine & cosine raised to a power */
		ParameterContainerT circ("circumferential");		
		{
			ParameterT m(ParameterT::Integer, "number_of_peaks");
			ParameterT k(ParameterT::Double, "concentration_k");
			ParameterT phi(ParameterT::Double, "location_phi");

			m.AddLimit(2, LimitT::LowerInclusive);
			k.AddLimit(-1.0e-13, LimitT::Lower);
			phi.AddLimit(-3.14159, LimitT::LowerInclusive);
			phi.AddLimit(3.14159, LimitT::UpperInclusive);
			
			circ.AddParameter(m);
			circ.AddParameter(k);
			circ.AddParameter(phi);
		}
		/* set the description */	
		circ.SetDescription("D(theta) = 1/(2 Pi I0(k)) exp(k*cos*(m*(theta-phi))) ");	
		choice->AddSub(circ);

		ParameterContainerT orthogonal("meridional");		
		{
			ParameterT m(ParameterT::Integer, "number_of_peaks");
			ParameterT k(ParameterT::Double, "concentration_k");
			ParameterT phi(ParameterT::Double, "location_phi");

			m.AddLimit(2, LimitT::LowerInclusive);
			k.AddLimit(-1.0e-13, LimitT::Lower);
			phi.AddLimit(-3.14159, LimitT::LowerInclusive);
			phi.AddLimit(3.14159, LimitT::UpperInclusive);
			
			orthogonal.AddParameter(m);
			orthogonal.AddParameter(k);
			orthogonal.AddParameter(phi);
		}
		orthogonal.SetDescription("D(theta) = 1/(2 Pi I0(k)) exp(k*cos*(m*(theta-phi))) ");	
		choice->AddSub(orthogonal);
		
		return(choice);
	}
}

/* accept parameter list */
void AnisoCornea::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FSFiberMatSplitT::TakeParameterList(list);
		
	const ParameterListT& matrix = list.GetListChoice(*this, "matrix_material_params");
	if (matrix.Name() == "Neo-Hookean")
	{
		fMu = matrix.GetParameter("shear_modulus");
		fKappa = matrix.GetParameter("bulk_modulus");
	}

	const ParameterListT& potential = list.GetListChoice(*this, "fibril_potential");
	if (potential.Name() == "fung_type")
	{
		double alpha = potential.GetParameter("alpha");
		double beta = potential.GetParameter("beta");
		fPotential = new FungType(alpha, beta);
		if (!fPotential) throw ExceptionT::kOutOfMemory;
	}
	
	const ParameterListT& distr = list.GetListChoice(*this, "fibril_distribution");
	if (distr.Name() == "circumferential")
	{
		fm = distr.GetParameter("number_of_peaks");
		fk = distr.GetParameter("concentration_k");
		fphi = distr.GetParameter("location_phi");
		 
		 fDType = kCircumferential;
	}
	else if (distr.Name() == "meridional")
	{
		fm = distr.GetParameter("number_of_peaks");
		fk = distr.GetParameter("concentration_k");
		fphi = distr.GetParameter("location_phi");
		 
		 fDType = kMeridional;
	}
		
	/* allocate memory */
	/*2D fiber stress and modulus*/
	fNumFibStress = dSymMatrixT::NumValues(fNumSD-1);
	fNumFibModuli = dSymMatrixT::NumValues(fNumFibStress);
	fFiberStretch.Dimension(fNumSD-1);
	fFiberStress.Dimension(fNumSD-1);
	fFiberMod.Dimension(fNumFibStress);

	/* point generator */
	int points = list.GetParameter("n_points");
	fCircle = new EvenSpacePtsT(points);
	Construct();
}

/***********************************************************************
 * Protected
 ***********************************************************************/

void AnisoCornea::ComputeDevMatrixStress(const dSymMatrixT& Cbar, dSymMatrixT& Sbar)
{
	/*2pdf{W}{Cbar_IJ} = mu*/
	
	Sbar[0] += fMu;
	Sbar[1] += fMu;
	Sbar[2] += fMu;

}

void AnisoCornea::ComputeDevMatrixMod(const dSymMatrixT& Cbar, dSymMatrixT& Stress, dMatrixT& Mod)
{
	/*2pdf{Sbar_IJ}{Cbar_KL} = mu*/
	
	Stress[0] += fMu;
	Stress[1] += fMu;
	Stress[2] += fMu;

	Mod += 0.0;

}

double AnisoCornea::ComputeVolMatrixStress(const double I3)
{
	/*p = 2pdf(Wvol)(I3)I3 = kappa/2(I3-1) C^-1*/
	return (0.5*fKappa*(I3-1));
}

double AnisoCornea::ComputeVolMatrixMod(const double I3)
{
	/*2pdf(p)(I3)*I3*/
	return(fKappa*I3);
}
	
/*computes integrated fiber stress in local frame*/
void AnisoCornea::ComputeFiberStress (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress)
{
	/*initialize pointers*/
	/* PK2 values in local frame formed by NT and IS orientations*/	
	FiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_12*/

	double *p1  = fStressTable(0);
	double *p2  = fStressTable(1);
	double *p3  = fStressTable(2);

	/* stretched bonds */
	ComputeLengths(FiberStretch);

	/* derivatives of the potential */
	fPotential->MapDFunction(fI4, fdU);

	/* initialize kernel pointers */
	double* pdU = fdU.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();
	
	/*integrate w.r.t in-plane orientation theta*/
	for (int i = 0; i < fI4.Length(); i++)
	{
		double factor = 2.0*(*pj++)*(*pdU++);
		s1 += factor*(*p1++);
		s2 += factor*(*p2++);
		s3 += factor*(*p3++);
	}

}
	
/*computes integrated moduli in local frame*/
void AnisoCornea::ComputeFiberMod (const dSymMatrixT& FiberStretch, dSymMatrixT& FiberStress, dMatrixT& FiberMod)
{
	/*initialize pointers*/
	/* PK2 and Modulus values in local coordinate frame formed by a1 (fNT) and a2 (fIS) */	
	FiberStress = 0.0;
	double& s1 = FiberStress[0]; /*sf_11*/
	double& s2 = FiberStress[1]; /*sf_22*/
	double& s3 = FiberStress[2]; /*sf_12*/
 
	FiberMod = 0.0;
	double& c11 = FiberMod[0]; /*cf_1111*/ 
	double& c22 = FiberMod[4]; /*cf_2222*/
	double& c33 = FiberMod[8]; /*cf_1212*/
	double& c23 = FiberMod[7]; /*cf_2212*/
	double& c13 = FiberMod[6]; /*cf_1112*/
	double& c12 = FiberMod[3]; /*cf_1122*/

	/*0  3  6
	  1  4  7
	  2  5  8*/

	/* stress */
	double* ps1  = fStressTable(0);
	double* ps2  = fStressTable(1);
	double* ps3  = fStressTable(2);

	/* modulus */
	double* pc11  = fModuliTable(0);
	double* pc22  = fModuliTable(1);
	double* pc33  = fModuliTable(2);
	double* pc23  = fModuliTable(3);
	double* pc13  = fModuliTable(4);
	double* pc12  = fModuliTable(5);

	/* stretched bonds */
	ComputeLengths(FiberStretch);
	
	/* derivatives of the potential */
	fPotential->MapDFunction(fI4, fdU);
	fPotential->MapDDFunction(fI4, fddU);	

	/* initialize kernel pointers */	
	double* pdU  = fdU.Pointer();
	double* pddU = fddU.Pointer();
	double* pj = fjacobians[fFSFiberMatSupport->CurrElementNumber()].Pointer();

		
	for (int i = 0; i < fI4.Length(); i++)
	{
		double sfactor =  2.0*(*pj)*(*pdU++);
		double cfactor = 4.0*(*pj++)*(*pddU++);

		s1 += sfactor*(*ps1++);
		s2 += sfactor*(*ps2++);
		s3 += sfactor*(*ps3++);

		c11 += cfactor*(*pc11++);
		c22 += cfactor*(*pc22++);
		c33 += cfactor*(*pc33++);

		c23 += cfactor*(*pc23++);
		c13 += cfactor*(*pc13++);
		c12 += cfactor*(*pc12++);
			
	}
		
	/*symmetric modulus*/
	FiberMod[1] = c12;
	FiberMod[2] = c13;
	FiberMod[5] = c23;
}

/* strained lengths in terms of the Lagrangian stretch eigenvalues */
void AnisoCornea::ComputeLengths(const dSymMatrixT& FiberStretch)
{	
	/*calculate fibril lengths                                                        *
	 *I4 = C*:M where M = cos^2 a1 x a1 + sin^2 a2 x a2 + sin cos (a1 x a2 + a2 x a1) */

	const double& C11 = FiberStretch[0];
	const double& C22 = FiberStretch[1];
	const double& C12 = FiberStretch[2];

	/* initialize kernel pointers */
	double* pl = fI4.Pointer();
	double* s1 = fStressTable(0);
	double* s2 = fStressTable(1);
	double* s3 = fStressTable(2);
		
	for (int i = 0; i < fI4.Length(); i++)
		*pl++ = C11*(*s1++) + C22*(*s2++) + 2.0*C12*(*s3++);
}

/***********************************************************************
* Private
***********************************************************************/
/* Initialize angle tables */
void AnisoCornea::Construct(void)
{
	/* fetch points */
	const dArray2DT& points = fCircle->CirclePoints(0.0);
	int numbonds = points.MajorDim();
		
	int nel = NumElements();     
	fjacobians.Dimension(nel);

	/* length table */
	fI4.Dimension(numbonds);

	/* potential tables */
	fU.Dimension(numbonds);
	fdU.Dimension(numbonds);
	fddU.Dimension(numbonds);

	/* STRESS angle tables - by associated stress component */
	fStressTable.Dimension(fNumFibStress, numbonds);
	  	
	/* MODULI angle tables - using Cauchy symmetry */ 	
	fModuliTable.Dimension(fNumFibModuli, numbonds);	

	/* set pointers */
	double *s1 = fStressTable(0);
	double *s2 = fStressTable(1);
	double *s3 = fStressTable(2);

	double *c11 = fModuliTable(0);
	double *c22 = fModuliTable(1);
	double *c33 = fModuliTable(2);
	double *c23 = fModuliTable(3);
	double *c13 = fModuliTable(4);
	double *c12 = fModuliTable(5);

	StringT name("distributions.dat");
	ofstreamT dist_out(name);


	/*modified von-mises distribution function [0,Pi/2]. integral[D(theta),{theta,-Pi,Pi}] = 2 Pi */
	if (fDType == kMeridional)
	{			 		
		/*integrate dtemp - min to normalize distribution function*/
		dArrayT& jac = fjacobians[0];
		jac.Dimension(numbonds);
		
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(fk);

		for (int i = 0; i < numbonds; i++)
		{
			/*function*/
			double theta = angles[i];
			double y = fk*cos(fm*(theta-fphi));

			double f = exp(y)/(I0*2.0*Pi);
			jac[i] = f * weights[i];
		}
		
		for (int el = 1; el < nel; el++)
		{
			fjacobians[el].Dimension(numbonds);			
			fjacobians[el] = jac;
		}
	}
	/*modified von-mises distribution function [0,Pi]. integral[D(theta),{theta,-Pi,Pi}] = 2 Pi */
	else if(fDType == kCircumferential)
	{
		dArrayT xc(3);
		const dArray2DT& coordinates = fFSFiberMatSupport->InitialCoordinates();

		fFSFiberMatSupport->TopElement();
		int ielm = fFSFiberMatSupport->CurrElementNumber();
	
		const dArrayT& angles = fCircle->CircleAngles(0.0);
		const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/
		const dArrayT& weights = fCircle->Jacobians();
		
		double I0 = besseli0(fk);

		while (fFSFiberMatSupport->NextElement()) 
		{
			int ielm = fFSFiberMatSupport->CurrElementNumber();
		
			/* calculate element centroid */
			iArrayT nodes = fFSFiberMatSupport->CurrentElement()->NodesX();
			int nen = NumElementNodes();       
			xc = 0.0;
			for (int i = 0; i < nen; i++) 
			{
				for (int j = 0; j < coordinates.MinorDim(); j++)
				xc[j] += coordinates(nodes[i],j);
			}
			xc /= nen;
		
			int inormal = 2; // HACK
			double x1 = xc[perm[inormal][1]];
			double x2 = xc[perm[inormal][2]];
		
			// projected polar coordinates Assume center is at (0,0) 
			double r = sqrt(x1*x1 + x2*x2);
			double xi = atan2(x2,x1);
			
		//shift not yet implemented//
			xi -= 0.5*Pi;
			xi -= fphi;
				
			dArrayT& jac = fjacobians[ielm];			
			jac.Dimension(numbonds);
			
			const dArrayT& angles = fCircle->CircleAngles(0.0);
			const dArray2DT& points = fCircle->CirclePoints(0.0);  /*set jacobians*/
			const dArrayT& weights = fCircle->Jacobians();
			
			double I0 = besseli0(fk);

			for (int i = 0; i < numbonds; i++)
			{
				/*function*/
				double theta = angles[i];
				double y = fk*cos(fm*(theta-xi));

				double f = exp(y)/(I0*2.0*Pi);
				jac[i] = f * weights[i];
			}
/*
			int elem = CurrElementNumber();
			if (elem == 1791 || elem == 2333 || elem == 2674 || elem == 2659 || elem == 2612 || elem == 1771)
			{
				cout << "\nelem: "<<CurrElementNumber();
				cout << "\ncentroid: "<<xc;
				cout << "\nangle: "<<xi;
				cout << "\nangle2: "<<xi;
				cout << "\njac: "<<jac;
				cout << "\nangle: "<<angles;
			}
*/
		}
	}	
	for (int i = 0; i < numbonds; i++)
	{
		/* direction cosines */
		const double *xsi = points(i);
		double cosi = xsi[0];
		double sini = xsi[1];
		
		/* stress angle tables */
		s1[i] = cosi*cosi;      
		s2[i] = sini*sini;
		s3[i] = sini*cosi;
	
		/* moduli angle tables */
		c11[i] = s1[i]*s1[i];
		c22[i] = s2[i]*s2[i];
		c33[i] = s3[i]*s3[i];
		c23[i] = s2[i]*s3[i];
		c13[i] = s1[i]*s3[i];
		c12[i] = s2[i]*s1[i];
	}
}
#endif /*VIB_MATERIAL*/

