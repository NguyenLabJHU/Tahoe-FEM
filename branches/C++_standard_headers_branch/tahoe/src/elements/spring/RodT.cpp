/* $Id: RodT.cpp,v 1.35 2005-11-06 00:37:58 paklein Exp $ */
/* created: paklein (10/22/1996) */
#include "RodT.h"

#include <math.h>
#include "ifstreamT.h"
#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"
#include "ParameterContainerT.h"

/* material types */
#include "LinearSpringT.h"
#include "LJSpringT.h"

using namespace Tahoe;

/* constructors */
RodT::RodT(const ElementSupportT& support):
	ElementBaseT(support),
	fCurrMaterial(NULL),
	fLocAcc(LocalArrayT::kAcc),
	fInstKE(0.0),
	fInstPE(0.0),
	fInstTotalE(0.0),
	fInstTemp(0.0),
	fInstPressure(0.0),
	fAvgKE(0.0),
	fAvgPE(0.0),
	fAvgTotalE(0.0),
	fAvgTemp(0.0),
	fAvgPressure(0.0),
	fSumKE(0.0),
	fSumPE(0.0),
	fSumTotalE(0.0),
	fSumTemp(0.0),
	fSumPressure(0.0),
	fLocVel(LocalArrayT::kVel),
	fOutputDiagnostic(false),
	fKb(1.38054)
{
	/* class name */
	SetName("spring_element");
}

/* form of tangent matrix */
GlobalT::SystemTypeT RodT::TangentType(void) const
{
	/* special case */
	if (fIntegrator->Order() > 0 &&
	    fIntegrator->ImplicitExplicit() ==  eIntegratorT::kExplicit)
		return GlobalT::kDiagonal;
	else
		return GlobalT::kSymmetric;
}

/* NOT implemented. Returns an zero force vector */
void RodT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
	const char caller[] = "RodT::AddNodalForce";

	/* different field */
	if (field.FieldName() != Field().FieldName()) return;

	/* set components and weights */
	double constKd = 0.0;
	int formKd = fIntegrator->FormKd(constKd);
	
	/* coordinates arrays */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	
	/* point forces */
	dArrayT f0(NumDOF(), fRHS.Pointer());
	dArrayT f1(NumDOF(), fRHS.Pointer() + NumDOF());

	Top();
	while (NextElement())
	{
		/* node numbers */
		const iArrayT& nodes = CurrentElement().NodesX();
		int n0 = nodes[0];
		int n1 = nodes[1];
		if (node == n0 || node == n1)
		{
			/* reference bond */
			fBond0.DiffOf(init_coords(n1), init_coords(n0));
			double l0 = fBond0.Magnitude();

			/* current bond */
			fBond.DiffOf(curr_coords(n1), curr_coords(n0));
			double l = fBond.Magnitude();
		
			/* bond force magnitude */
			double dU = fCurrMaterial->DPotential(l, l0);
			double f_by_l = 0.0;
			if (fabs(l) > kSmall)
				f_by_l = constKd*dU/l;
			else if (fabs(dU) > kSmall)
				ExceptionT::GeneralFail(caller, "bond %d has length but %g force", CurrElementNumber(), dU);

			/* particle forces (extra -1 since moved to the RHS) */
			f0.SetToScaled(f_by_l, fBond);
			
			/* assemble */
			if (node == n0)
				force -= f0;
			else
				force += f0;
		}
	}
}

/* returns the energy as defined by the derived class types */
double RodT::InternalEnergy(void)
{
	/* coordinates arrays */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();

	double energy = 0.0;
	Top();
	while (NextElement())
	{
		/* node numbers */
		const iArrayT& nodes = CurrentElement().NodesX();
		int n0 = nodes[0];
		int n1 = nodes[1];
	
		/* reference bond */
		fBond0.DiffOf(init_coords(n1), init_coords(n0));

		/* current bond */
		fBond.DiffOf(curr_coords(n1), curr_coords(n0));
		
		/* form element stiffness */
		energy += fCurrMaterial->Potential(fBond.Magnitude(), fBond0.Magnitude());
	}
	return energy;
}

/* writing output */
void RodT::RegisterOutput(void)
{
	/* block ID's */
	ArrayT<StringT> block_ID(fBlockData.Length());
	for (int i = 0; i < block_ID.Length(); i++)
		block_ID[i] = fBlockData[i].ID();

	/* set output specifier */
	ArrayT<StringT> e_labels;
	OutputSetT output_set(GeometryT::kLine, block_ID, fConnectivities, 
		Field().Labels(), e_labels, ChangingGeometry());
		
	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

void RodT::WriteOutput(void)
{
	/* get list of nodes used by the group */
	iArrayT nodes_used;
	NodesUsed(nodes_used);

	/* temp space for group displacements */
	dArray2DT disp(nodes_used.Length(), NumDOF());
	
	/* collect group displacements */
	disp.RowCollect(nodes_used, Field()[0]);

	/* send */
	dArray2DT e_values;
	ElementSupport().WriteOutput(fOutputID, disp, e_values);
}

/* compute specified output parameter and send for smoothing */
void RodT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

void RodT::CloseStep(void)
{
	/* inherited */
	ElementBaseT::CloseStep();

	if (fOutputDiagnostic) {
		Top();
		const FieldT& field = Field();
		while (NextElement()) {
			ComputeHardyStress();
			if (field.Order() > 0)
			{
				/* get velocities */
				SetLocalU(fLocVel);
	
				/* compute MD quantities of interest */
				ComputeInstPE();
				ComputeInstKE();
				ComputeInstTotalE();
				ComputeInstTemperature();
				//ComputeInstPressure();
      		}
  		}
  		ComputeAvgPE();
  		ComputeAvgKE();
  		ComputeAvgTotalE();
  		ComputeAvgTemperature();
  		//ComputeAvgPressure()
  		PrintMDToFile();
  	}
}

/* describe the parameters needed by the interface */
void RodT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ElementBaseT::DefineParameters(list);

	/* output diagnostic data */
	ParameterT output_diagnostic(fOutputDiagnostic, "output_diagnostic");
	output_diagnostic.SetDefault(fOutputDiagnostic);
	list.AddParameter(output_diagnostic);
}

/* information about subordinate parameter lists */
void RodT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ElementBaseT::DefineSubs(sub_list);

	/* element block/constitutive specification */
	sub_list.AddSub("spring_element_block", ParameterListT::OnePlus);
}

/* return the description of the given inline subordinate parameter list */
ParameterInterfaceT* RodT::NewSub(const StringT& name) const
{
	if (name == "spring_element_block")
	{
		ParameterContainerT* block = new ParameterContainerT(name);
		
		/* list of element block ID's (defined by ElementBaseT) */
		block->AddSub("block_ID_list", ParameterListT::Once);
	
		/* choice of materials lists (inline) */
		block->AddSub("spring_material_choice", ParameterListT::Once, true);
	
		/* set this as source of subs */
		block->SetSubSource(this);
		
		return block;
	} 
	else if (name == "spring_material_choice") 
	{
		ParameterContainerT* mat_choice = new ParameterContainerT(name);
		mat_choice->SetListOrder(ParameterListT::Choice);
	
		/* mass parameter */
		ParameterT mass(ParameterT::Double, "mass");
		mass.AddLimit(0, LimitT::LowerInclusive);

		/* linear spring */
		ParameterContainerT linear("spring_linear");
		linear.AddParameter(mass);
		ParameterT k(ParameterT::Double, "stiffness");
		k.AddLimit(0, LimitT::LowerInclusive);
		linear.AddParameter(k);
		mat_choice->AddSub(linear);

		/* Lennard-Jones spring */
		ParameterContainerT LJ("spring_LJ");
		LJ.AddParameter(mass);
		ParameterT eps(ParameterT::Double, "epsilon");
		eps.AddLimit(0, LimitT::LowerInclusive);
		LJ.AddParameter(eps);
		ParameterT sigma(ParameterT::Double, "sigma");
		sigma.AddLimit(0, LimitT::LowerInclusive);
		LJ.AddParameter(sigma);
		mat_choice->AddSub(LJ);
	
		return mat_choice;
	}
	else /* inherited */
		return ElementBaseT::NewSub(name);
}

/* accept parameter list */
void RodT::TakeParameterList(const ParameterListT& list)
{
	const char caller[] = "RodT::TakeParameterList";

	/* inherited */
	ElementBaseT::TakeParameterList(list);

	/* output diagnostic data */
	fOutputDiagnostic = list.GetParameter("output_diagnostic");

	/* dimensions */
	int nsd = NumSD();
	int ndof = NumDOF();
	
	/* constant matrix needed to calculate stiffness */
	fOneOne.Dimension(fLHS);
	dMatrixT one(ndof);
	one.Identity();
	fOneOne.SetBlock(0, 0, one);
	fOneOne.SetBlock(ndof, ndof, one);
	one *= -1;
	fOneOne.SetBlock(0, ndof, one);
	fOneOne.SetBlock(ndof, 0, one);
	
	/* work space */
	fBond.Dimension(nsd);
	fBond0.Dimension(nsd);
	fHardyStress.Dimension(nsd);
	fHardyHeatFlux.Dimension(nsd);	

	/* echo material properties */
//	ReadMaterialData(ElementSupport().Input());	
//	WriteMaterialData(ElementSupport().Output());

	/* get form of tangent */
	GlobalT::SystemTypeT type = TangentType();
	
	/* set form of element stiffness matrix */
	if (type == GlobalT::kSymmetric)
		fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
	else if (type == GlobalT::kNonSymmetric)
		fLHS.SetFormat(ElementMatrixT::kNonSymmetric);
	else if (type == GlobalT::kDiagonal)
		fLHS.SetFormat(ElementMatrixT::kDiagonal);

	/* local arrays */
	const FieldT& field = Field();
	fLocVel.Dimension(NumElementNodes(), ndof);
	if (field.Order() > 0) Field().RegisterLocal(fLocVel);
	fLocAcc.Dimension(NumElementNodes(), ndof);
	if (fIntegrator->Order() == 2) Field().RegisterLocal(fLocAcc);
	fNEE_vec.Dimension(fLocAcc.Length());

	/* determine the nodes used strictly based on those in the connectivities */
	NodesUsed(fGroupNodes);
	
	/* number of element blocks */
	int num_blocks = fConnectivities.Length();
	fMaterialsList.Dimension(num_blocks);
	for (int i = 0; i < num_blocks; i++)
	{
		/* block parameters */
		const ParameterListT* block_info = list.FindList("_element_block", i);
		if (!block_info) ExceptionT::GeneralFail(caller, "could not resolve instance %d of block information", i+1);

		/* material parameters */
		const ParameterListT* mat_info = block_info->ListChoice(*this, "spring_material_choice");
		if (!mat_info) ExceptionT::GeneralFail(caller, "could not resolve \"spring_material_choice\"", i+1);
		if (mat_info->Name() == "spring_linear")
		{
			/* parameters */
			double mass = mat_info->GetParameter("mass");
			double k = mat_info->GetParameter("stiffness");

			/* construct */
			fMaterialsList[i] = new LinearSpringT(mass, k);
		}
		else if (mat_info->Name() == "spring_LJ")
		{
			/* parameters */
			double mass = mat_info->GetParameter("mass");
			double eps = mat_info->GetParameter("epsilon");
			double sigma = mat_info->GetParameter("sigma");

			/* construct */
			fMaterialsList[i] = new LJSpringT(mass, eps, sigma);
		}
		else
			ExceptionT::GeneralFail(caller, "could not resolve type \"%s\"",
				mat_info->Name().Pointer());
	}
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct the element stiffness matrix */
void RodT::LHSDriver(GlobalT::SystemTypeT)
{
	/* time integration dependent */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fIntegrator->FormK(constK);
	int formM = fIntegrator->FormM(constM);

	/* coordinates arrays */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	
	/* use RHS as temp space */
	dArrayT v0(NumDOF(), fRHS.Pointer());
	dArrayT v1(NumDOF(), fRHS.Pointer() + NumDOF());

	Top();
	while (NextElement())
	{
		/* particle mass */
		double mass = fCurrMaterial->Mass();

		/* node numbers */
		const iArrayT& nodes = CurrentElement().NodesX();
		int n0 = nodes[0];
		int n1 = nodes[1];

		/* form element stiffness */
		if (formK) {		

			/* reference bond */
			fBond0.DiffOf(init_coords(n1), init_coords(n0));
			double l0 = fBond0.Magnitude();

			/* current bond */
			v1.DiffOf(curr_coords(n1), curr_coords(n0));
			double l = v1.Magnitude();
			v0.SetToScaled(-1.0, v1);
			
			/* bond force and stiffness */
			double dU = fCurrMaterial->DPotential(l, l0);
			double ddU = fCurrMaterial->DDPotential(l, l0);

			/* allow zero length springs */
			if (fabs(l) > kSmall) 
			{
				/* 1st term */
				fLHS.Outer(fRHS, fRHS);
				fLHS *= constK*(ddU - dU/l)/(l*l);
		
				/* 2nd term */
				fLHS.AddScaled(constK*dU/l, fOneOne);
			}
			else /* limit as l->0 */
				fLHS.AddScaled(constK*ddU, fOneOne);
		} 
		else 
			fLHS = 0.0;
	
		/* mass contribution */
		if (formM) fLHS.PlusIdentity(constM*mass);
	
		/* add to global equations */
		AssembleLHS();
	}
}

/* construct the element force vectors */
void RodT::RHSDriver(void)
{
	const char caller[] = "RodT::RHSDriver";

	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fIntegrator->FormMa(constMa);
	int formKd = fIntegrator->FormKd(constKd);
	
	/* coordinates arrays */
	const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
	const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
	
	/* point forces */
	dArrayT f0(NumDOF(), fRHS.Pointer());
	dArrayT f1(NumDOF(), fRHS.Pointer() + NumDOF());

	Top();
	while (NextElement())
	{
		/* node numbers */
		const iArrayT& nodes = CurrentElement().NodesX();
		int n0 = nodes[0];
		int n1 = nodes[1];
	
		/* reference bond */
		fBond0.DiffOf(init_coords(n1), init_coords(n0));
		double l0 = fBond0.Magnitude();

		/* current bond */
		fBond.DiffOf(curr_coords(n1), curr_coords(n0));
		double l = fBond.Magnitude();
		
		/* bond force magnitude */
		double dU = fCurrMaterial->DPotential(l, l0);
		double f_by_l = 0.0;
		if (fabs(l) > kSmall)
			f_by_l = constKd*dU/l;
		else if (fabs(dU) > kSmall)
			ExceptionT::GeneralFail(caller, "bond %d has length but %g force", CurrElementNumber(), dU);

		/* particle forces (extra -1 since moved to the RHS) */
		f0.SetToScaled(f_by_l, fBond);
		f1.SetToScaled(-1.0, f0);

		/* inertial force */
		if (formMa) {		
			SetLocalU(fLocAcc);
			fLocAcc.ReturnTranspose(fNEE_vec);
			fRHS.AddScaled(-constMa*fCurrMaterial->Mass(), fNEE_vec);
		}

		/* add to global equations */
		AssembleRHS();
	}
}

/* load next element */
bool RodT::NextElement(void)
{
	bool result = ElementBaseT::NextElement();
	
	/* initialize element calculation */
	if (result)
		fCurrMaterial = fMaterialsList[CurrentElement().MaterialNumber()];
	
	return result;
}

/***********************************************************************
* Private
***********************************************************************/

/* Below are functions implementing Hardy ideas deriving continuum 
 * measures from MD notions */

void RodT::ComputeHardyStress(void)
{
  double constKd = 0.0;
	
  /* components dicated by the algorithm */
  int formKd = fIntegrator->FormKd(constKd);
  
  /* coordinates arrays */
  const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
  const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
  
  /* point forces */
  dArrayT f0(NumDOF(), fRHS.Pointer());
  dArrayT f1(NumDOF(), fRHS.Pointer() + NumDOF());
  
  /* node numbers */
  const iArrayT& nodes = CurrentElement().NodesX();
  int n0 = nodes[0];
  int n1 = nodes[1];
	
  /* reference bond */
  fBond0.DiffOf(init_coords(n1), init_coords(n0));
  double l0 = fBond0.Magnitude();

  /* current bond */
  fBond.DiffOf(curr_coords(n1), curr_coords(n0));
  double l = fBond.Magnitude();
		
  /* bond force magnitude */
  double dU = fCurrMaterial->DPotential(l, l0);
	
  /* particle forces */
  f0.SetToScaled(constKd*dU/l, fBond);

  /* Potential part of Hardy stress */
  fHardyStress.Outer(fBond,f0);
  fHardyStress *= -.5;

  /* Kinetic part of Hardy stress */
  const FieldT& field = Field();
  if (field.Order() > 0) {
   
    dArrayT vel;
    dMatrixT kinstress(NumSD());
    const dArray2DT& velocities = field[1];
    for (int i = 0; i < fGroupNodes.Length(); i++)
      {
	velocities.RowAlias(fGroupNodes[i], vel);
	kinstress.Outer(vel,vel);
	kinstress *= .5;
	kinstress *= fCurrMaterial->Mass();
	fHardyStress += kinstress;
      }
  }
}

void RodT::ComputeHardyHeatFlux(void) { }

/* Below are MD related computational functions */

void RodT::ComputeInstKE(void)
{
	/* computes the instantaneous kinetic energy of the system of atoms */
	double& ke = fInstKE;
	double& total = fSumKE;
  
	ke = 0.0;
	const FieldT& field = Field();
	if (field.Order() > 0) {
   
		dArrayT vel;
		const dArray2DT& velocities = field[1];
		for (int i = 0; i < fGroupNodes.Length(); i++)
  		{
  			velocities.RowAlias(fGroupNodes[i], vel);
  			ke += dArrayT::Dot(vel,vel);
		}
		ke *= fCurrMaterial->Mass()/2.0;
	}
	total += ke;
}

void RodT::ComputeAvgKE(void)
{
	double& tempavg = fAvgKE;
	int nstep = ElementSupport().StepNumber();
	if (nstep) /* only for positive steps */
  		tempavg = fSumKE / nstep;
	else
		tempavg = 0.0;
}

void RodT::ComputeInstPE(void)
{
  /* computes the potential energy of the system of atoms */
  
  /* coordinates arrays */
  const dArray2DT& init_coords = ElementSupport().InitialCoordinates();
  const dArray2DT& curr_coords = ElementSupport().CurrentCoordinates();
  double& pe = fInstPE;
  double& total = fSumPE;
  pe = 0.0;

  /* node numbers */
  const iArrayT& nodes = CurrentElement().NodesX();
  int n0 = nodes[0];
  int n1 = nodes[1];
  
  /* reference bond */
  fBond0.DiffOf(init_coords(n1), init_coords(n0));
  
  /* current bond */
  fBond.DiffOf(curr_coords(n1), curr_coords(n0));

  pe += fCurrMaterial->Potential(fBond.Magnitude(), fBond0.Magnitude());
  total += pe;
}

void RodT::ComputeAvgPE(void)
{
	double& tempavg = fAvgPE;
	int nstep = ElementSupport().StepNumber();	
  	if (nstep > 0) /* only for positive steps */
  		tempavg = fSumPE / nstep;
	else
		tempavg = 0.0;
}

void RodT::ComputeInstTotalE(void)
{
  /* computes instantaneous total energy = kinetic energy + potential energy of the system */
  double& totale = fInstTotalE;
  double& total = fSumTotalE;
  totale = 0.0;
  totale = fInstKE + fInstPE;
  total += totale;
}

void RodT::ComputeAvgTotalE(void)
{
  	double& tempavg = fAvgTotalE;
	int nstep = ElementSupport().StepNumber();  	
	if (nstep > 0) /* only for positive steps */
		tempavg = fSumTotalE / nstep;
	else
		tempavg = 0.0;
}

void RodT::ComputeInstTemperature(void)
{
  /* computes instantaneous temperature of the atomic system */
  double& temp = fInstTemp;
  temp = 0.0;
  temp = 2 * fInstKE / (fKb * NumSD() * fGroupNodes.Length());
}

void RodT::ComputeAvgTemperature(void)
{
  double& temp = fAvgTemp;
  temp = 2 * fAvgKE / (fKb * NumSD() * fGroupNodes.Length());
}

void RodT::ComputeInstPressure(void)
{
  /* computes the instantaneous pressure of the atomic system */

}

void RodT::ComputeAvgPressure(void)
{
  /* computes the average pressure of the atomic system */

}

int RodT::PrintMDToFile(void)
{
	/* print MD quantities (temperature/energy/pressure) to an output file */
	ofstreamT out;
	int d_width = OutputWidth(out, &fAvgTotalE); 
	if (ElementSupport().StepNumber() == 1)
  	{
  		out.open("MD.out");
		if (!out.is_open()) {
			cout << "Cannot open MD.out file.\n";
			return 1;
		}

  		out << setw(d_width) << "Timestep"
		    << setw(d_width) << "InstKE"  
  		    << setw(d_width) << "InstPE" 
  		    << setw(d_width) << "InstTemp" 
  		    << setw(d_width) << "InstTotalE" 
  		    << setw(d_width) << "AvgKE" 
  		    << setw(d_width) << "AvgPE" 
  		    << setw(d_width) << "AvgTemp" 
  		    << setw(d_width) << "AvgTotalE" << endl;
	}
	else /* appending to existing file */
	{
  		out.open_append("MD.out");
		if (!out.is_open()) {
			cout << "Cannot open MD.out file.\n";
			return 1;
		}
	
		out << setw(d_width) << ElementSupport().StepNumber()
		    << setw(d_width) << fInstKE 
		    << setw(d_width) << fInstPE 
		    << setw(d_width) << fInstTemp 
		    << setw(d_width) << fInstTotalE 
		    << setw(d_width) << fAvgKE 
		    << setw(d_width) << fAvgPE 
		    << setw(d_width) << fAvgTemp 
		    << setw(d_width) << fAvgTotalE 
		    << setw(d_width) << endl;
	}
	out.close();
	return 0;
}
