/* $Id: ParticlePairT.cpp,v 1.3 2002-11-21 01:11:14 paklein Exp $ */
/* created: paklein (10/22/1996) */

#include "ParticleT.h"

#include <math.h>

#include "fstreamT.h"
#include "eControllerT.h"
#include "OutputSetT.h"
#include "dArray2DT.h"

/* interaction types */
#include "LennardJones612.h"
#include "ParabolaT.h"
#include "SmithFerrante.h"

/* constructors */
ParticleT::ParticleT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field)
{
	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);
}

/* initialization */
void ParticleT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();
	
	/* constant matrix needed to calculate stiffness */
#if 0
	fOneOne.Dimension(fLHS);
	dMatrixT one(NumDOF());
	one.Identity();
	fOneOne.SetBlock(0, 0, one);
	fOneOne.SetBlock(NumDOF(), NumDOF(), one);
	one *= -1;
	fOneOne.SetBlock(0, NumDOF(), one);
	fOneOne.SetBlock(NumDOF(), 0, one);

	/* bond vector */
	fBond.Dimension(NumSD());
#endif
	
	/* echo material properties */
	ReadMaterialData(ElementSupport().Input());	
	WriteMaterialData(ElementSupport().Output());
}

/* form of tangent matrix */
GlobalT::SystemTypeT ParticleT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* NOT implemented. Returns an zero force vector */
void ParticleT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

void ParticleT::WriteOutput(IOBaseT::OutputModeT mode)
{
//TEMP - not handling general output modes yet
	if (mode != IOBaseT::kAtInc)
	{
		cout << "\n ContinuumElementT::WriteOutput: only handling \"at increment\"\n"
		     <<   "     print mode. SKIPPING." << endl;
		return;
	}

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
void ParticleT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, do nothing
}

/***********************************************************************
* Protected
***********************************************************************/

/* construct the element stiffness matrix */
void ParticleT::LHSDriver(void)
{
	/* time integration dependent */
	double constK = 0.0;
	double constM = 0.0;
	int formK = fController->FormK(constK);
	int formM = fController->FormM(constM);
	
	/* particle mass */
	double mass = 0.0; //fCurrMaterial->Mass();
	
	Top();
	while ( NextElement() )
	{
		/* initialize */
		fLHS = 0.0;
		
		/* local arrays */
		SetLocalX(fLocInitCoords);
		SetLocalU(fLocDisp);
		
		/* form element stiffness */
		if (formK) ElementStiffness(constK);
	
		/* mass contribution */
		if (formM) fLHS.PlusIdentity(mass);
	
		/* add to global equations */
		AssembleLHS();
	}
}

/* construct the element force vectors */
void ParticleT::RHSDriver(void)
{
	/* set components and weights */
	double constMa = 0.0;
	double constKd = 0.0;
	
	/* components dicated by the algorithm */
	int formMa = fController->FormMa(constMa);
	int formKd = fController->FormKd(constKd);
	
//TEMP - inertia term in residual
if (formMa) {
	cout << "\n ParticleT::RHSDriver: M*a term not implemented" << endl;
	throw eGeneralFail;
}

	/* run through pairs */
	Top();
	while (NextElement()) /* would be faster not to use Top-Next */
	{
		/* local displacement */
		SetLocalU(fLocDisp);
	
		if (fLocDisp.AbsMax() > 0.0 || fCurrMaterial->HasInternalStrain())
		{
			/* initialize */
			fRHS = 0.0;
	
			/* local coordinates */
			SetLocalX(fLocInitCoords);

			/* form element force */
			ElementForce(-1.0);
	
			/* add to global equations */
			AssembleRHS();
		}
	}
}

/* load next element */
bool ParticleT::NextElement(void)
{
	bool result = ElementBaseT::NextElement();
	
	/* initialize element calculation */
	if (result)
		fCurrMaterial = fMaterialsList[CurrentElement().MaterialNumber()];
	
	return result;
}
	
/* element data */
void ParticleT::ReadMaterialData(ifstreamT& in)
{
	/* allocate space */
	int	nummaterials;
	in >> nummaterials;
	fInteractions.Dimension(nummaterials);
	fPointMass.Dimension(nummaterials);

	/* read data */
	for (int i = 0; i < nummaterials; i++)
	{
		int matnum, matcode;
		in >> matnum;
		matnum--;

		/* associated mass */
		double mass = = -1;
		in >> mass;
		if (mass < 0) throw eBadInputValue;
		fPointMass[matnum] = mass;
		
		/* add to the list of materials */
		in >> matcode;
		switch (matcode)
		{
			case kLennardJones:			
			{
				double a = -1;
				in >> a;
				if (a < 0) throw eBadInputValue;
				fInteractions[matnum] = new LennardJones612(a);
				break;
			}	
			case kSmithFerrante:
			{
				double A = -1, B = -1, L = -1;
				in >> A >> B >> L;
				if (A < 0 || B < 0 || L < 0) throw eBadInputValue;
				fInteractions[matnum] = new SmithFerrante(A, B, L);
				break;
			}
			case kQuadratic:
			{
				double k = -1;
				in >> k;
				if (k < 0) throw eBadInputValue;
				fInteractions[matnum] = new ParabolaT(k);
				break;
			}
			default:
				cout << "\n ParticleT::ReadMaterialData: unknown material type\n" << endl;
				throw eBadInputValue;
		}
	}
}

void ParticleT::WriteMaterialData(ostream& out) const
{
	out << "\n Particle Set Data:\n";
	out << " Number of particle types. . . . . . . . . . = " << fInteractions.Length() << '\n';
	for (int i = 0; i < fInteractions.Length(); i++)
	{
		out << "\n Type number . . . . . . . . . . . . . . . . = " << i+1 << '\n';
		out << " Mass. . . . . . . . . . . . . . . . . . . . = " << fPointMass[i] << '\n';
		out << " Pair interaction:\n";
		fInteractions[i]->PrintName(out);
		fInteractions[i]->Print(out);
	}
}
