/* $Id: SWDiamondT.cpp,v 1.12 2004-06-17 07:41:37 paklein Exp $ */
/* created: paklein (03/19/1997) */

#include "SWDiamondT.h"

#include <math.h>
#include <iomanip.h>

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "FindNeighbor23T.h"
#include "OutputSetT.h"

/* element group parameters */

using namespace Tahoe;

const int knsd  = 3;
const int kSWMaxNeighbors0 = 4; //max neighbors in undeformed state

/* constructor */
SWDiamondT::SWDiamondT(const ElementSupportT& support, const FieldT& field):
	ElementBaseT(support, field),
	fK_3Body(fLHS),
	fF_3Body(fRHS),
	fK_2Body(ElementMatrixT::kSymmetricUpper),
	fLocX_3Body(LocalArrayT::kInitCoords, 3, knsd),
	fLocd_3Body(LocalArrayT::kDisp, 3, NumDOF()),
	List_3Body(fElementCards),
	fLocX_2Body(LocalArrayT::kInitCoords, 2, knsd),
	fLocd_2Body(LocalArrayT::kDisp, 2, NumDOF()),
	fHessian_3Body(3)
{
	/* check base class initializations */
	if (NumSD() != knsd) throw ExceptionT::kGeneralFail;

	/* set matrix format */
	fLHS.SetFormat(ElementMatrixT::kSymmetricUpper);

	/* set base class data */
//	fNumElemEqnos = 3*NumDOF();

	/* allocate memory */
	fK_3Body.Dimension(3*NumDOF());
	fF_3Body.Dimension(3*NumDOF());

	fK_2Body.Dimension(2*NumDOF());
	fF_2Body.Dimension(2*NumDOF());

	/* register local arrays */
	ElementSupport().RegisterCoordinates(fLocX_3Body);
	ElementSupport().RegisterCoordinates(fLocX_2Body);
	
	Field().RegisterLocal(fLocd_3Body);
	Field().RegisterLocal(fLocd_2Body);
	
	//TEMP
	ReadMaterialData(ElementSupport().Input());	
}

/* initialization */
void SWDiamondT::Initialize(void)
{
	/* inherited */
	ElementBaseT::Initialize();

	/* material parameters */
	//ReadMaterialData(ElementSupport().Input());	
	WriteMaterialData(ElementSupport().Output());
}

/* form of tangent matrix */
GlobalT::SystemTypeT SWDiamondT::TangentType(void) const
{
	return GlobalT::kSymmetric;
}

/* NOT implemented. Returns an zero force vector */
void SWDiamondT::AddNodalForce(const FieldT& field, int node, dArrayT& force)
{
#pragma unused(field)
#pragma unused(node)
#pragma unused(force)
}

/* returns the energy as defined by the derived class types */
double SWDiamondT::InternalEnergy(void)
{
	double energy = 0.0;

	/* 3 body contribution */
	List_3Body.Top();
	while ( Next3Body() )
	{
		/* local arrays */
		SetLocalX(fLocX_3Body);
		SetLocalU(fLocd_3Body);
		
		/* form element stiffness */
		energy += Energy3Body();
	}

	/* 2 body contribution */
	List_2Body.Top();
	while ( Next2Body() )
	{
		/* local arrays */
		SetLocalX(fLocX_2Body);
		SetLocalU(fLocd_2Body);
		
		/* form element stiffness */
		energy += Energy2Body();
	}

	return(energy);
}

/* append element equations numbers to the list */
void SWDiamondT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
	AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_2)

	/* set local equations numbers */
	Field().SetLocalEqnos(fNodes_3Body, fEqnos_3Body);
	Field().SetLocalEqnos(fNodes_2Body, fEqnos_2Body);

	/* add to list */
	eq_1.Append(&fEqnos_3Body); // 3-body interactions include 2-body
}
	
/* writing output */
void SWDiamondT::RegisterOutput(void)
{
	/* set labels */
	ArrayT<StringT> n_labels(3), e_labels;
	n_labels[0] = "D_X";
	n_labels[1] = "D_Y";
	n_labels[2] = "D_Z";		

	/* set block IDs, since fBlockData is nil */
	ArrayT<StringT> block_ID (1);
	block_ID[0] = "1";

	/* set output specifier */
	ArrayT<const iArray2DT*> output_connects_list(1);
	output_connects_list[0] = &fOutputConnects;
	OutputSetT output_set(GeometryT::kPoint, block_ID, output_connects_list, n_labels, e_labels, false);

	/* register and get output ID */
	fOutputID = ElementSupport().RegisterOutput(output_set);
}

void SWDiamondT::WriteOutput(void)
{
	/* calculate output values */
	dArray2DT n_values(fNodesUsed.Length(), NumDOF()), e_values;
	n_values.RowCollect(fNodesUsed, Field()[0]); /* displacements */

	/* send to output */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

/* compute specified output parameter and send for smoothing */
void SWDiamondT::SendOutput(int kincode)
{
#pragma unused(kincode)
	//TEMP: for now, does nothing
}
/***********************************************************************
* Protected
***********************************************************************/

/* construct the element stiffness matrix */
void SWDiamondT::LHSDriver(GlobalT::SystemTypeT)
{
	/* 3 body contribution */
	List_3Body.Top();
	while (List_3Body.Next()) //assume ALL will contribute
	{
		/* initialize */
		fK_3Body = 0.0;

		/* current interaction */
		ElementCardT& card = List_3Body.Current();
	
		/* local arrays */
		fLocX_3Body.SetLocal(card.NodesX());
		fLocd_3Body.SetLocal(card.NodesU());
						
		/* form element stiffness */
		Stiffness3Body();
	
		/* add to global equations */
		ElementSupport().AssembleLHS(Group(), fK_3Body, card.Equations());
	}

	/* 2 body contribution */
	List_2Body.Top();
	while (List_2Body.Next()) //assume ALL will contribute
	{
		/* initialize */
		fK_2Body = 0.0;

		/* current interaction */
		ElementCardT& card = List_2Body.Current();

		/* local arrays */
		fLocX_2Body.SetLocal(card.NodesX());
		fLocd_2Body.SetLocal(card.NodesX());
		
		/* form element stiffness */
		Stiffness2Body();
	
		/* add to global equations */
		ElementSupport().AssembleLHS(Group(), fK_2Body, card.Equations());
	}
}

/* construct the element force vectors */
void SWDiamondT::RHSDriver(void)
{
	/* 3 body contribution */
	List_3Body.Top();
	while (List_3Body.Next()) //assume ALL will contribute
	{
		/* initialize */
		fF_3Body = 0.0;

		/* current interaction */
		ElementCardT& card = List_3Body.Current();
	
		/* local arrays */
		fLocX_3Body.SetLocal(card.NodesX());
		fLocd_3Body.SetLocal(card.NodesU());
						
		/* form element stiffness */
		Force3Body();
	
		/* add to global equations */
		ElementSupport().AssembleRHS(Group(), fF_3Body, card.Equations());
	}

	/* 2 body contribution */
	List_2Body.Top();
	while (List_2Body.Next()) //assume ALL will contribute
	{
		/* initialize */
		fF_2Body = 0.0;

		/* current interaction */
		ElementCardT& card = List_2Body.Current();

		/* local arrays */
		fLocX_2Body.SetLocal(card.NodesX());
		fLocd_2Body.SetLocal(card.NodesX());
		
		/* form element stiffness */
		Force2Body();
	
		/* add to global equations */
		ElementSupport().AssembleRHS(Group(), fF_2Body, card.Equations());
	}
}

/* print element group data */
void SWDiamondT::PrintControlData(ostream& out) const
{
#pragma unused(out)
}
	
/* element data */
void SWDiamondT::ReadMaterialData(ifstreamT& in)
{
	/* unit scaling */
	in >> feps;	if (feps <= 0.0) throw ExceptionT::kBadInputValue;

	/* 2 body potential */
	in >> fA;		if (fA     <= 0.0) throw ExceptionT::kBadInputValue;
	in >> fdelta;	if (fdelta <= 0.0) throw ExceptionT::kBadInputValue;
	
	/* 3 body potential */
	in >> fgamma;	if (fgamma  <= 0.0) throw ExceptionT::kBadInputValue;
	in >> flambda;	if (flambda <= 0.0) throw ExceptionT::kBadInputValue;
	
	in >> frcut;	if (frcut <= 0.0) throw ExceptionT::kBadInputValue;		
	in >> fa;		if (fa    <= 0.0) throw ExceptionT::kBadInputValue;

	/* set B factor */
	double a0 = pow(2.0,1.0/6.0);
	fB =-(fdelta*pow(a0,5))/(-(a0*fdelta) - 4.0*a0*a0 + 8.0*a0*frcut - 4.0*frcut*frcut);
}

void SWDiamondT::WriteMaterialData(ostream& out) const
{
	out << "\n Material Set Data:\n";

	/* unit scaling */
	out << "    epsilon = " << feps << '\n';		

	/* 2 body potential */
	out << " 2 body terms:\n";
	out << "          A = " << fA << '\n';
	out << "          B = " << fB << "   **COMPUTED\n";
	out << "      delta = " << fdelta << '\n';
	
	/* 3 body potential */
	out << " 3 body terms:\n";
	out << "      gamma = " << fgamma << '\n';
	out << "     lambda = " << flambda << '\n';
	
	out << " Lattice scaling and cut-off terms:\n";
	out << "      r_cut = " << frcut << '\n';
	out << "         a0 = " << fa << '\n';
}

void SWDiamondT::EchoConnectivityData(ifstreamT& in, ostream& out)
{
	int num_nodes_used;
	in >> num_nodes_used;
	if (num_nodes_used != -1 &&
	    num_nodes_used  < 1) throw ExceptionT::kBadInputValue;

	/* neighbor distance */
	double tolerance = 1.01*pow(2.0,1.0/6.0)*fa;

	/* read nodes used */
	if (num_nodes_used == -1) //use ALL nodes
	{
		/* current coordinates */
		const dArray2DT& coords = ElementSupport().CurrentCoordinates();
	
		/* connector */
		FindNeighbor23T Connector(coords, kSWMaxNeighbors0);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(fNodes_2Body, fNodes_3Body, tolerance);
		
		/* set nodes used */
		fNodesUsed.Dimension(coords.MajorDim());
		fNodesUsed.SetValueToPosition();
	}
	else                      //only use specified nodes
	{
		/* read specified nodes */
		fNodesUsed.Dimension(num_nodes_used);
		in >> fNodesUsed;

		/* echo data */
		out << "\n Nodes used : \n\n";
		int* pnodes    = fNodesUsed.Pointer();
		int  linecount = 0;
		for (int i = 0; i < num_nodes_used; i++)
		{
			out << setw(kIntWidth) << *pnodes++;
			
			if (++linecount == 5)
			{
				out << '\n';
				linecount = 0;
			}
			else
				out << "   ";
		}
		if (linecount != 0) out << '\n';

		/* connector */
		FindNeighbor23T Connector(fNodesUsed, ElementSupport().CurrentCoordinates(),
									kSWMaxNeighbors0);
	
		/* connect nodes - dimensions lists */
		Connector.GetNeighors(fNodes_2Body, fNodes_3Body, tolerance);		
	}
	
	/* connectivities of "point elements" */
	fOutputConnects.Set(fNodesUsed.Length(), 1, fNodesUsed.Pointer());
	
	/* set element equation and node lists */
	ConfigureElementData();

	/* print connectivity data */
	PrintConnectivityData(out);
}

/* call AFTER 2 and 3 body node lists are set */
void SWDiamondT::ConfigureElementData(void)
{
	int num3 = fNodes_3Body.MajorDim();
	int num2 = fNodes_2Body.MajorDim();

	int neq3body = 3*NumDOF();
	int neq2body = 2*NumDOF();

	/* allocate memory */
	List_3Body.Dimension(num3);
	List_2Body.Dimension(num2);

	fEqnos_3Body.Dimension(num3, neq3body);
	fEqnos_2Body.Dimension(num2, neq2body);

	/* set 3 body element data */
	for (int i = 0; i < num3; i++)	
	{
		/* node and equation numbers */			
		(List_3Body[i].NodesX()).Set(3, fNodes_3Body(i));		
		(List_3Body[i].Equations()).Set(neq3body, fEqnos_3Body(i));
	}
	
	/* set 2 body element data */
	for (int j = 0; j < num2; j++)	
	{
		/* node and equation numbers */			
		(List_2Body[j].NodesX()).Set(2, fNodes_2Body(j));		
		(List_2Body[j].Equations()).Set(neq2body, fEqnos_2Body(j));
	}
	
	/* set base class connectivity and equations data */
	fConnectivities.Dimension(1);
	fConnectivities[0] = &fNodes_3Body;
	fEqnos.Dimension(1);
	fEqnos[0].Alias(fEqnos_3Body);
}

/* element list increment */
bool SWDiamondT::Next2Body(void) { return List_2Body.Next(); }
bool SWDiamondT::Next3Body(void) { return List_3Body.Next(); }

/* element calculations */
double SWDiamondT::Energy3Body(void)
{	
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;
	double z12, z13, z14, z15, z16, z17, z18;
	
	z1 = fLocd_3Body(0,0);
	z2 = fLocd_3Body(0,1);
	z3 = fLocd_3Body(0,2);
	z4 = fLocd_3Body(1,0);
	z5 = fLocd_3Body(1,1);
	z6 = fLocd_3Body(1,2);
	z7 = fLocd_3Body(2,0);
	z8 = fLocd_3Body(2,1);
	z9 = fLocd_3Body(2,2);
	z10 = fLocX_3Body(0,0);
	z11 = fLocX_3Body(0,1);
	z12 = fLocX_3Body(0,2);
	z13 = fLocX_3Body(1,0);
	z14 = fLocX_3Body(1,1);
	z15 = fLocX_3Body(1,2);
	z16 = fLocX_3Body(2,0);
	z17 = fLocX_3Body(2,1);
	z18 = fLocX_3Body(2,2);
	z13 = -z13;
	z14 = -z14;
	z15 = -z15;
	z4 = -z4;
	z5 = -z5;
	z6 = -z6;
	z1 = z1 + z10 + z13 + z4;
	z2 = z11 + z14 + z2 + z5;
	z3 = z12 + z15 + z3 + z6;
	z4 = z13 + z16 + z4 + z7;
	z5 = z14 + z17 + z5 + z8;
	z6 = z15 + z18 + z6 + z9;
	z7 = z1*z1;
	z8 = pow(z2,2.);
	z9 = pow(z3,2.);
	z1 = z1*z4;
	z4 = pow(z4,2.);
	z2 = z2*z5;
	z5 = pow(z5,2.);
	z3 = z3*z6;
	z6 = pow(z6,2.);
	z7 = z7 + z8 + z9;
	z1 = z1 + z2 + z3;
	z2 = z4 + z5 + z6;
	z3 = pow(z7,-0.5);
	z4 = pow(z7,0.5);
	z5 = pow(z2,-0.5);
	z2 = pow(z2,0.5);
	z1 = z1*z3*z5;

	//U2[z4, z2, z1]

	return( U3body(z4, z2, z1) );
}

void SWDiamondT::Force3Body(void)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22;
	double z23, z24, z25, z26, z27, z28, z29, z30, z31, z32, z33, z34;
	double z35, z36;

	z1 = fLocd_3Body(0,0);
	z2 = fLocd_3Body(0,1);
	z3 = fLocd_3Body(0,2);
	z4 = fLocd_3Body(1,0);
	z5 = fLocd_3Body(1,1);
	z6 = fLocd_3Body(1,2);
	z7 = fLocd_3Body(2,0);
	z8 = fLocd_3Body(2,1);
	z9 = fLocd_3Body(2,2);
	z10 = fLocX_3Body(0,0);
	z11 = fLocX_3Body(0,1);
	z12 = fLocX_3Body(0,2);
	z13 = fLocX_3Body(1,0);
	z14 = fLocX_3Body(1,1);
	z15 = fLocX_3Body(1,2);
	z16 = fLocX_3Body(2,0);
	z17 = fLocX_3Body(2,1);
	z18 = fLocX_3Body(2,2);
	z19 = -1.*z1;
	z20 = -1.*z10;
	z21 = -1.*z11;
	z22 = -1.*z12;
	z23 = -1.*z13;
	z13 = 2.*z13;
	z24 = -1.*z14;
	z14 = 2.*z14;
	z25 = -1.*z15;
	z15 = 2.*z15;
	z26 = -1.*z16;
	z27 = -1.*z17;
	z28 = -1.*z18;
	z29 = -1.*z2;
	z30 = -1.*z3;
	z31 = -1.*z4;
	z4 = 2.*z4;
	z32 = -1.*z5;
	z5 = 2.*z5;
	z33 = -1.*z6;
	z6 = 2.*z6;
	z34 = -1.*z7;
	z35 = -1.*z8;
	z36 = -1.*z9;
	z1 = z1 + z10 + z23 + z31;
	z2 = z11 + z2 + z24 + z32;
	z3 = z12 + z25 + z3 + z33;
	z4 = z13 + z19 + z20 + z26 + z34 + z4;
	z5 = z14 + z21 + z27 + z29 + z35 + z5;
	z6 = z15 + z22 + z28 + z30 + z36 + z6;
	z7 = z16 + z23 + z31 + z7;
	z8 = z17 + z24 + z32 + z8;
	z9 = z18 + z25 + z33 + z9;
	z10 = pow(z1,2.);
	z11 = pow(z2,2.);
	z12 = pow(z3,2.);
	z13 = z1*z7;
	z14 = pow(z7,2.);
	z15 = z2*z8;
	z16 = pow(z8,2.);
	z17 = z3*z9;
	z18 = pow(z9,2.);
	z10 = z10 + z11 + z12;
	z11 = z13 + z15 + z17;
	z12 = z14 + z16 + z18;
	z13 = pow(z10,-1.5);
	z14 = pow(z10,-0.5);
	z10 = pow(z10,0.5);
	z15 = pow(z12,-1.5);
	z16 = pow(z12,-0.5);
	z12 = pow(z12,0.5);
	z17 = -1.*z11;
	z18 = z11*z14*z15;
	z19 = z1*z11*z13*z16;
	z20 = z11*z13*z16*z2;
	z21 = z11*z13*z16*z3;
	z22 = z1*z14*z16;
	z23 = z14*z16*z2;
	z24 = z14*z16*z3;
	z4 = z14*z16*z4;
	z5 = z14*z16*z5;
	z6 = z14*z16*z6;
	z25 = z14*z16*z7;
	z26 = z14*z16*z8;
	z27 = z14*z16*z9;
	z11 = z11*z14*z16;
	z15 = z14*z15*z17;
	z28 = z1*z13*z16*z17;
	z29 = z13*z16*z17*z2;
	z13 = z13*z16*z17*z3;
	z17 = z18*z7;
	z30 = z18*z8;
	z18 = z18*z9;
	z31 = z15*z7;
	z32 = z15*z8;
	z15 = z15*z9;
	z25 = z25 + z28;
	z26 = z26 + z29;
	z13 = z13 + z27;
	z4 = z17 + z19 + z4;
	z5 = z20 + z30 + z5;
	z6 = z18 + z21 + z6;
	z17 = z22 + z31;
	z18 = z23 + z32;
	z15 = z15 + z24;
	z19 = Dc12U3body(z10,z12,z11);
	z20 = Dr2U3body(z10,z12,z11);
	z10 = Dr1U3body(z10,z12,z11);
	z11 = -1.*z19;
	z6 = z11*z6;
	z12 = z11*z17;
	z17 = z11*z18;
	z15 = z11*z15;
	z18 = -1.*z1*z10*z14;
	z1 = z1*z10*z14;
	z19 = -1.*z10*z14*z2;
	z2 = z10*z14*z2;
	z21 = -1.*z10*z14*z3;
	z3 = z10*z14*z3;
	z10 = -1.*z16*z20*z7;
	z7 = z16*z20*z7;
	z14 = -1.*z16*z20*z8;
	z8 = z16*z20*z8;
	z22 = -1.*z16*z20*z9;
	z9 = z16*z20*z9;
	z16 = z11*z25;
	z20 = z11*z26;
	z13 = z11*z13;
	z4 = z11*z4;
	z5 = z11*z5;
	z10 = z10 + z12;
	z11 = z14 + z17;
	z12 = z15 + z22;
	z3 = z3 + z6 + z9;
	z6 = z16 + z18;
	z9 = z19 + z20;
	z13 = z13 + z21;
	z1 = z1 + z4 + z7;
	z2 = z2 + z5 + z8;
	
	//List(z6,z9,z13,z1,z2,z3,z10,z11,z12);
	
	fF_3Body[0] = z6;
	fF_3Body[1] = z9;
	fF_3Body[2] = z13;
	fF_3Body[3] = z1;
	fF_3Body[4] = z2;
	fF_3Body[5] = z3;
	fF_3Body[6] = z10;
	fF_3Body[7] = z11;
	fF_3Body[8] = z12;
}

void SWDiamondT::Stiffness3Body(void)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27, z28, z29, z30, z31, z32, z33, z34, z35, z36;
	double z37, z38, z39, z40, z41, z42, z43, z44, z45, z46, z47, z48;
	double z49, z50, z51, z52, z53, z54, z55, z56, z57, z58, z59, z60;
	double z61, z62, z63, z64, z65, z66, z67, z68, z69, z70, z71, z72;
	double z73, z74, z75, z76, z77, z78, z79, z80, z81, z82, z83, z84;
	double z85, z86, z87, z88, z89, z90, z91, z92, z93, z94, z95, z96;
	double z97, z98, z99, z100, z101, z102, z103, z104, z105, z106, z107, z108;
	double z109, z110, z111, z112, z113, z114, z115, z116, z117, z118, z119, z120;
	double z121, z122, z123, z124, z125, z126, z127, z128, z129, z130, z131, z132;
	double z133, z134, z135, z136, z137, z138, z139, z140, z141, z142, z143, z144;
	double z145, z146, z147, z148, z149, z150, z151, z152, z153, z154, z155, z156;
	double z157, z158, z159, z160, z161, z162, z163, z164, z165, z166, z167, z168;
	double z169, z170, z171, z172, z173, z174, z175, z176, z177, z178, z179, z180;
	double z181, z182, z183, z184, z185, z186, z187, z188;

	z1 = fLocd_3Body(0,0);
	z2 = fLocd_3Body(0,1);
	z3 = fLocd_3Body(0,2);
	z4 = fLocd_3Body(1,0);
	z5 = fLocd_3Body(1,1);
	z6 = fLocd_3Body(1,2);
	z7 = fLocd_3Body(2,0);
	z8 = fLocd_3Body(2,1);
	z9 = fLocd_3Body(2,2);
	z10 = fLocX_3Body(0,0);
	z11 = fLocX_3Body(0,1);
	z12 = fLocX_3Body(0,2);
	z13 = fLocX_3Body(1,0);
	z14 = fLocX_3Body(1,1);
	z15 = fLocX_3Body(1,2);
	z16 = fLocX_3Body(2,0);
	z17 = fLocX_3Body(2,1);
	z18 = fLocX_3Body(2,2);
	z19 = -1.*z1;
	z20 = -1.*z10;
	z21 = -1.*z11;
	z22 = -1.*z12;
	z23 = -1.*z13;
	z13 = 2.*z13;
	z24 = -1.*z14;
	z14 = 2.*z14;
	z25 = -1.*z15;
	z15 = 2.*z15;
	z26 = -1.*z16;
	z27 = -1.*z17;
	z28 = -1.*z18;
	z29 = -1.*z2;
	z30 = -1.*z3;
	z31 = -1.*z4;
	z4 = 2.*z4;
	z32 = -1.*z5;
	z5 = 2.*z5;
	z33 = -1.*z6;
	z6 = 2.*z6;
	z34 = -1.*z7;
	z35 = -1.*z8;
	z36 = -1.*z9;
	z1 = z1 + z10 + z23 + z31;
	z2 = z11 + z2 + z24 + z32;
	z3 = z12 + z25 + z3 + z33;
	z4 = z13 + z19 + z20 + z26 + z34 + z4;
	z5 = z14 + z21 + z27 + z29 + z35 + z5;
	z6 = z15 + z22 + z28 + z30 + z36 + z6;
	z7 = z16 + z23 + z31 + z7;
	z8 = z17 + z24 + z32 + z8;
	z9 = z18 + z25 + z33 + z9;
	z10 = pow(z1,2.);
	z11 = pow(z2,2.);
	z12 = pow(z3,2.);
	z13 = z1*z7;
	z14 = pow(z7,2.);
	z15 = z2*z8;
	z16 = pow(z8,2.);
	z17 = z3*z9;
	z18 = pow(z9,2.);
	z19 = z10 + z11 + z12;
	z20 = z13 + z15 + z17;
	z21 = z14 + z16 + z18;
	z22 = pow(z19,-2.5);
	z23 = pow(z19,-1.5);
	z24 = pow(z19,-0.5);
	z19 = pow(z19,0.5);
	z25 = pow(z21,-2.5);
	z26 = pow(z21,-1.5);
	z27 = pow(z21,-0.5);
	z21 = pow(z21,0.5);
	z28 = z20*z7;
	z29 = -3.*z20*z24*z25;
	z30 = z20*z24*z25;
	z31 = -1.*z1*z20*z23*z26*z8;
	z32 = z1*z20*z23*z26*z8;
	z33 = -1.*z20*z23*z26*z3*z8;
	z34 = z20*z23*z26*z3*z8;
	z35 = -1.*z1*z20*z23*z26*z9;
	z36 = z1*z20*z23*z26*z9;
	z37 = -1.*z2*z20*z23*z26*z9;
	z38 = z2*z20*z23*z26*z9;
	z39 = -1.*z13*z20*z23*z26;
	z40 = z13*z20*z23*z26;
	z41 = -1.*z15*z20*z23*z26;
	z42 = z15*z20*z23*z26;
	z43 = -1.*z17*z20*z23*z26;
	z44 = z17*z20*z23*z26;
	z45 = -1.*z2*z24*z26*z7;
	z46 = z2*z24*z26*z7;
	z47 = -1.*z24*z26*z3*z7;
	z48 = z24*z26*z3*z7;
	z49 = -1.*z24*z26*z4*z7;
	z50 = 2.*z24*z26*z4*z7;
	z51 = -1.*z24*z26*z5*z7;
	z52 = z24*z26*z5*z7;
	z53 = -1.*z24*z26*z6*z7;
	z54 = z24*z26*z6*z7;
	z55 = -1.*z1*z24*z26*z8;
	z56 = z1*z24*z26*z8;
	z57 = -1.*z24*z26*z3*z8;
	z58 = z24*z26*z3*z8;
	z59 = -1.*z24*z26*z4*z8;
	z60 = z24*z26*z4*z8;
	z61 = -1.*z24*z26*z5*z8;
	z62 = 2.*z24*z26*z5*z8;
	z63 = -1.*z24*z26*z6*z8;
	z64 = z24*z26*z6*z8;
	z65 = -1.*z24*z26*z7*z8;
	z66 = z24*z26*z7*z8;
	z67 = -1.*z1*z24*z26*z9;
	z68 = z1*z24*z26*z9;
	z69 = -1.*z2*z24*z26*z9;
	z70 = z2*z24*z26*z9;
	z71 = -1.*z24*z26*z4*z9;
	z72 = z24*z26*z4*z9;
	z73 = -1.*z24*z26*z5*z9;
	z74 = z24*z26*z5*z9;
	z75 = -1.*z24*z26*z6*z9;
	z76 = 2.*z24*z26*z6*z9;
	z77 = -1.*z24*z26*z7*z9;
	z78 = z24*z26*z7*z9;
	z79 = -1.*z24*z26*z8*z9;
	z80 = z24*z26*z8*z9;
	z81 = -2.*z13*z24*z26;
	z82 = z13*z24*z26;
	z83 = -1.*z14*z24*z26;
	z84 = z14*z24*z26;
	z85 = -2.*z15*z24*z26;
	z86 = z15*z24*z26;
	z87 = -1.*z16*z24*z26;
	z88 = z16*z24*z26;
	z89 = -2.*z17*z24*z26;
	z90 = z17*z24*z26;
	z91 = -1.*z18*z24*z26;
	z92 = z18*z24*z26;
	z93 = -1.*z20*z24*z26;
	z94 = z8*z93;
	z95 = z9*z93;
	z96 = z20*z24*z26;
	z97 = z8*z96;
	z98 = z9*z96;
	z99 = -3.*z1*z2*z20*z22*z27;
	z100 = 3.*z1*z2*z20*z22*z27;
	z101 = -3.*z1*z20*z22*z27*z3;
	z102 = 3.*z1*z20*z22*z27*z3;
	z103 = -3.*z2*z20*z22*z27*z3;
	z104 = 3.*z2*z20*z22*z27*z3;
	z105 = -3.*z10*z20*z22*z27;
	z106 = 3.*z10*z20*z22*z27;
	z107 = -3.*z11*z20*z22*z27;
	z108 = 3.*z11*z20*z22*z27;
	z109 = -3.*z12*z20*z22*z27;
	z22 = 3.*z12*z20*z22*z27;
	z110 = -1.*z1*z2*z23*z27;
	z111 = z1*z2*z23*z27;
	z112 = -1.*z1*z23*z27*z3;
	z113 = z1*z23*z27*z3;
	z114 = -1.*z2*z23*z27*z3;
	z115 = z2*z23*z27*z3;
	z116 = -1.*z1*z23*z27*z4;
	z117 = 2.*z1*z23*z27*z4;
	z118 = -1.*z2*z23*z27*z4;
	z119 = z2*z23*z27*z4;
	z120 = -1.*z23*z27*z3*z4;
	z121 = z23*z27*z3*z4;
	z122 = -1.*z1*z23*z27*z5;
	z123 = z1*z23*z27*z5;
	z124 = -1.*z2*z23*z27*z5;
	z125 = 2.*z2*z23*z27*z5;
	z126 = -1.*z23*z27*z3*z5;
	z127 = z23*z27*z3*z5;
	z128 = -1.*z1*z23*z27*z6;
	z129 = z1*z23*z27*z6;
	z130 = -1.*z2*z23*z27*z6;
	z131 = z2*z23*z27*z6;
	z132 = -1.*z23*z27*z3*z6;
	z133 = 2.*z23*z27*z3*z6;
	z134 = -1.*z2*z23*z27*z7;
	z135 = z2*z23*z27*z7;
	z136 = -1.*z23*z27*z3*z7;
	z137 = z23*z27*z3*z7;
	z138 = -1.*z1*z23*z27*z8;
	z139 = z1*z23*z27*z8;
	z140 = -1.*z23*z27*z3*z8;
	z141 = z23*z27*z3*z8;
	z142 = -1.*z1*z23*z27*z9;
	z143 = z1*z23*z27*z9;
	z144 = -1.*z2*z23*z27*z9;
	z145 = z2*z23*z27*z9;
	z146 = -1.*z10*z23*z27;
	z147 = z10*z23*z27;
	z148 = -1.*z11*z23*z27;
	z149 = z11*z23*z27;
	z150 = -1.*z12*z23*z27;
	z151 = z12*z23*z27;
	z152 = -2.*z13*z23*z27;
	z13 = z13*z23*z27;
	z153 = -2.*z15*z23*z27;
	z15 = z15*z23*z27;
	z154 = -2.*z17*z23*z27;
	z17 = z17*z23*z27;
	z155 = -1.*z20*z23*z27;
	z156 = z1*z155;
	z157 = z155*z2;
	z158 = z155*z3;
	z159 = z20*z23*z27;
	z160 = z1*z159;
	z161 = z159*z2;
	z162 = z159*z3;
	z163 = -1.*z24*z27;
	z164 = z24*z27;
	z165 = 2.*z164;
	z166 = z1*z164;
	z167 = z164*z2;
	z168 = z164*z3;
	z4 = z164*z4;
	z5 = z164*z5;
	z6 = z164*z6;
	z169 = z164*z7;
	z170 = z164*z8;
	z171 = z164*z9;
	z20 = z164*z20;
	z172 = -3.*z24*z25*z28;
	z173 = z172*z8;
	z172 = z172*z9;
	z25 = z24*z25*z28;
	z174 = z25*z8;
	z174 = 3.*z174;
	z25 = 3.*z25*z9;
	z175 = -1.*z2*z23*z26*z28;
	z176 = z2*z23*z26*z28;
	z177 = -1.*z23*z26*z28*z3;
	z178 = z23*z26*z28*z3;
	z179 = -1.*z24*z26*z28;
	z28 = z24*z26*z28;
	z180 = z29*z8;
	z180 = z180*z9;
	z181 = z14*z29;
	z182 = z16*z29;
	z29 = z18*z29;
	z183 = z30*z8;
	z183 = 3.*z183*z9;
	z184 = 3.*z14*z30;
	z185 = 3.*z16*z30;
	z30 = 3.*z18*z30;
	z186 = 2.*z40;
	z187 = 2.*z42;
	z188 = 2.*z44;
	z134 = z100 + z134 + z138;
	z136 = z102 + z136 + z142;
	z138 = z104 + z140 + z144;
	z140 = z106 + z152 + z155;
	z142 = z108 + z153 + z155;
	z144 = z154 + z155 + z22;
	z94 = z167 + z94;
	z95 = z168 + z95;
	z5 = z161 + z5 + z97;
	z6 = z162 + z6 + z98;
	z97 = z156 + z169;
	z98 = z157 + z170;
	z152 = z158 + z171;
	z118 = z118 + z139 + z175 + z66 + z99;
	z139 = z110 + z176 + z65;
	z120 = z101 + z120 + z143 + z177 + z78;
	z143 = z112 + z178 + z77;
	z153 = z166 + z179;
	z4 = z160 + z28 + z4;
	z28 = z183 + z57 + z69;
	z57 = z184 + z81 + z93;
	z69 = z185 + z85 + z93;
	z81 = z30 + z89 + z93;
	z62 = z108 + z125 + z155 + z165 + z185 + z187 + z62 + z93;
	z22 = z133 + z155 + z165 + z188 + z22 + z30 + z76 + z93;
	z30 = z122 + z135 + z31 + z66 + z99;
	z65 = z110 + z32 + z65;
	z66 = z103 + z126 + z145 + z33 + z80;
	z33 = z115 + z180 + z33 + z63 + z70;
	z63 = z114 + z34 + z79;
	z70 = z101 + z128 + z137 + z35 + z78;
	z76 = z112 + z36 + z77;
	z77 = z103 + z130 + z141 + z37 + z80;
	z37 = z115 + z180 + z37 + z58 + z73;
	z58 = z114 + z38 + z79;
	z34 = z104 + z127 + z131 + z183 + z34 + z38 + z64 + z74;
	z13 = z105 + z116 + z13 + z159 + z163 + z39 + z84;
	z38 = z146 + z164 + z40 + z83;
	z15 = z107 + z124 + z15 + z159 + z163 + z41 + z88;
	z40 = z149 + z163 + z182 + z41 + z61 + z86 + z96;
	z41 = z148 + z164 + z42 + z87;
	z17 = z109 + z132 + z159 + z163 + z17 + z43 + z92;
	z29 = z151 + z163 + z29 + z43 + z75 + z90 + z96;
	z42 = z150 + z164 + z44 + z91;
	z31 = z111 + z173 + z31 + z46 + z59;
	z43 = z25 + z47 + z67;
	z35 = z113 + z172 + z35 + z48 + z71;
	z39 = z147 + z163 + z181 + z39 + z49 + z82 + z96;
	z44 = z106 + z117 + z155 + z165 + z184 + z186 + z50 + z93;
	z46 = z111 + z173 + z175 + z51 + z56;
	z32 = z100 + z119 + z123 + z174 + z176 + z32 + z52 + z60;
	z47 = z113 + z172 + z177 + z53 + z68;
	z25 = z102 + z121 + z129 + z178 + z25 + z36 + z54 + z72;
	z36 = z174 + z45 + z55;
	
//	z45 = Dc12U3body(z19,z21,z20);
//	z48 = DDc12U3body(z19,z21,z20);
//	z49 = Dr2U3body(z19,z21,z20);
//	z50 = Dr2Dc12U3body(z19,z21,z20);
//	z51 = DDr2U3body(z19,z21,z20);
//	z52 = Dr1U3body(z19,z21,z20);
//	z53 = Dr1Dc12U3body(z19,z21,z20);
//	z54 = Dr1Dr2U3body(z19,z21,z20);
//	z19 = DDr1U3body(z19,z21,z20);

	//pak (05/22/97)
	//Z'd version of the Hessian of the 3-body potential
	//doesn't really seem to speed things up.

	/* 1st derivatives */
	z45 = Dc12U3body(z19,z21,z20);
	z49 = Dr2U3body(z19,z21,z20);
	z52 = Dr1U3body(z19,z21,z20);

	/* compute Hessian components */
	ComputeHessian_3Body(fHessian_3Body, z19, z21, z20);
	
	/* assign values */
	z48 = fHessian_3Body(2,2); 		//DDc12U3body(z19,z21,z20);
	z50 = fHessian_3Body(1,2);		//Dr2Dc12U3body(z19,z21,z20);
	z51 = fHessian_3Body(1,1);		//DDr2U3body(z19,z21,z20);
	z53 = fHessian_3Body(0,2);		//Dr1Dc12U3body(z19,z21,z20);
	z54 = fHessian_3Body(0,1);		//Dr1Dr2U3body(z19,z21,z20);
	z19 = fHessian_3Body(0,0);		//DDr1U3body(z19,z21,z20);
	
	z20 = z134*z45;
	z21 = z136*z45;
	z55 = z138*z45;
	z56 = z140*z45;
	z59 = z142*z45;
	z60 = z144*z45;
	z61 = z118*z45;
	z64 = z139*z45;
	z67 = z120*z45;
	z68 = z143*z45;
	z28 = z28*z45;
	z57 = z45*z57;
	z69 = z45*z69;
	z71 = z45*z81;
	z62 = z45*z62;
	z22 = z22*z45;
	z30 = z30*z45;
	z65 = z45*z65;
	z66 = z45*z66;
	z33 = z33*z45;
	z63 = z45*z63;
	z70 = z45*z70;
	z72 = z45*z76;
	z73 = z45*z77;
	z37 = z37*z45;
	z58 = z45*z58;
	z34 = z34*z45;
	z13 = z13*z45;
	z38 = z38*z45;
	z15 = z15*z45;
	z40 = z40*z45;
	z41 = z41*z45;
	z17 = z17*z45;
	z29 = z29*z45;
	z42 = z42*z45;
	z31 = z31*z45;
	z43 = z43*z45;
	z35 = z35*z45;
	z39 = z39*z45;
	z44 = z44*z45;
	z46 = z45*z46;
	z32 = z32*z45;
	z47 = z45*z47;
	z25 = z25*z45;
	z36 = z36*z45;
	z45 = z48*z94;
	z74 = z48*z95;
	z75 = z48*z5;
	z76 = z48*z6;
	z77 = z48*z97;
	z78 = z48*z98;
	z79 = z152*z48;
	z80 = z153*z48;
	z48 = z4*z48;
	z81 = z50*z94;
	z82 = z50*z95;
	z83 = z5*z50;
	z84 = z50*z6;
	z85 = z153*z50;
	z86 = z4*z50;
	z87 = z53*z94;
	z88 = z53*z95;
	z89 = z5*z53;
	z90 = z53*z6;
	z91 = z53*z97;
	z92 = z53*z98;
	z93 = z152*z53;
	z96 = z153*z53;
	z99 = z4*z53;
	z100 = -1.*z1*z2*z23*z52;
	z101 = z1*z2*z23*z52;
	z102 = -1.*z1*z23*z3*z52;
	z103 = z1*z23*z3*z52;
	z104 = -1.*z2*z23*z3*z52;
	z105 = z2*z23*z3*z52;
	z106 = -1.*z10*z23*z52;
	z10 = z10*z23*z52;
	z107 = -1.*z11*z23*z52;
	z11 = z11*z23*z52;
	z108 = -1.*z12*z23*z52;
	z12 = z12*z23*z52;
	z23 = -1.*z24*z52;
	z52 = z24*z52;
	z109 = -1.*z1*z24*z53;
	z110 = z1*z24*z53;
	z111 = -1.*z1*z24*z54;
	z112 = -1.*z1*z19*z24;
	z113 = z1*z19*z24;
	z114 = -1.*z2*z24*z53;
	z115 = z2*z24*z53;
	z116 = -1.*z2*z24*z54;
	z117 = -1.*z19*z2*z24;
	z118 = z19*z2*z24;
	z119 = -1.*z24*z3*z53;
	z53 = z24*z3*z53;
	z120 = -1.*z24*z3*z54;
	z121 = -1.*z19*z24*z3;
	z19 = z19*z24*z3;
	z122 = -1.*z26*z49*z7*z8;
	z123 = z26*z49*z7*z8;
	z124 = -1.*z26*z49*z7*z9;
	z125 = z26*z49*z7*z9;
	z126 = -1.*z26*z49*z8*z9;
	z127 = z26*z49*z8*z9;
	z128 = -1.*z14*z26*z49;
	z14 = z14*z26*z49;
	z129 = -1.*z16*z26*z49;
	z16 = z16*z26*z49;
	z130 = -1.*z18*z26*z49;
	z18 = z18*z26*z49;
	z26 = -1.*z27*z49;
	z49 = z27*z49;
	z131 = -1.*z27*z50*z7;
	z132 = z27*z50*z7;
	z133 = -1.*z27*z51*z7;
	z134 = z27*z51*z7;
	z135 = -1.*z27*z54*z7;
	z136 = z27*z54*z7;
	z137 = -1.*z27*z50*z8;
	z138 = z27*z50*z8;
	z139 = -1.*z27*z51*z8;
	z140 = z27*z51*z8;
	z141 = -1.*z27*z54*z8;
	z142 = z27*z54*z8;
	z143 = -1.*z27*z50*z9;
	z50 = z27*z50*z9;
	z144 = -1.*z27*z51*z9;
	z51 = z27*z51*z9;
	z145 = -1.*z27*z54*z9;
	z54 = z27*z54*z9;
	z77 = z110 + z77;
	z91 = z113 + z91;
	z78 = z115 + z78;
	z92 = z118 + z92;
	z53 = z53 + z79;
	z19 = z19 + z93;
	z48 = z109 + z131 + z48;
	z79 = z132 + z80;
	z80 = z111 + z133 + z86;
	z85 = z134 + z85;
	z86 = z112 + z135 + z99;
	z93 = z136 + z96;
	z75 = z114 + z137 + z75;
	z45 = z138 + z45;
	z83 = z116 + z139 + z83;
	z81 = z140 + z81;
	z89 = z117 + z141 + z89;
	z87 = z142 + z87;
	z76 = z119 + z143 + z76;
	z50 = z50 + z74;
	z74 = z120 + z144 + z84;
	z51 = z51 + z82;
	z82 = z121 + z145 + z90;
	z54 = z54 + z88;
	z77 = z77*z97;
	z84 = z78*z97;
	z78 = z78*z98;
	z88 = z53*z97;
	z90 = z53*z98;
	z53 = z152*z53;
	z96 = z48*z97;
	z99 = z48*z98;
	z109 = z152*z48;
	z48 = z4*z48;
	z110 = z5*z79;
	z111 = z6*z79;
	z112 = z79*z97;
	z113 = z79*z98;
	z114 = z152*z79;
	z115 = z153*z79;
	z79 = z4*z79;
	z116 = z5*z75;
	z117 = z75*z97;
	z118 = z75*z98;
	z119 = z152*z75;
	z75 = z4*z75;
	z120 = z45*z94;
	z121 = z45*z5;
	z131 = z45*z6;
	z132 = z45*z97;
	z133 = z45*z98;
	z134 = z152*z45;
	z135 = z153*z45;
	z45 = z4*z45;
	z136 = z5*z76;
	z137 = z6*z76;
	z138 = z76*z97;
	z139 = z76*z98;
	z140 = z152*z76;
	z76 = z4*z76;
	z94 = z50*z94;
	z95 = z50*z95;
	z5 = z5*z50;
	z6 = z50*z6;
	z97 = z50*z97;
	z98 = z50*z98;
	z141 = z152*z50;
	z142 = z153*z50;
	z4 = z4*z50;
	z50 = z1*z24*z91;
	z91 = z1*z24*z92;
	z92 = z2*z24*z92;
	z143 = z1*z19*z24;
	z144 = z19*z2*z24;
	z19 = z19*z24*z3;
	z145 = -1.*z1*z24*z86;
	z146 = z1*z24*z86;
	z147 = z2*z24*z86;
	z86 = z24*z3*z86;
	z148 = -1.*z1*z24*z93;
	z149 = z1*z24*z93;
	z150 = -1.*z2*z24*z93;
	z151 = z2*z24*z93;
	z152 = -1.*z24*z3*z93;
	z93 = z24*z3*z93;
	z153 = -1.*z1*z24*z89;
	z154 = z1*z24*z89;
	z155 = -1.*z2*z24*z89;
	z156 = z2*z24*z89;
	z89 = z24*z3*z89;
	z157 = -1.*z1*z24*z87;
	z158 = z1*z24*z87;
	z159 = -1.*z2*z24*z87;
	z160 = z2*z24*z87;
	z161 = -1.*z24*z3*z87;
	z87 = z24*z3*z87;
	z162 = -1.*z1*z24*z82;
	z163 = z1*z24*z82;
	z164 = -1.*z2*z24*z82;
	z165 = z2*z24*z82;
	z166 = -1.*z24*z3*z82;
	z82 = z24*z3*z82;
	z167 = -1.*z1*z24*z54;
	z1 = z1*z24*z54;
	z168 = -1.*z2*z24*z54;
	z2 = z2*z24*z54;
	z169 = -1.*z24*z3*z54;
	z3 = z24*z3*z54;
	z24 = -1.*z27*z7*z80;
	z54 = -1.*z27*z7*z85;
	z80 = z27*z7*z85;
	z170 = -1.*z27*z7*z83;
	z171 = -1.*z27*z7*z81;
	z172 = z27*z7*z81;
	z173 = -1.*z27*z7*z74;
	z174 = -1.*z27*z51*z7;
	z7 = z27*z51*z7;
	z175 = -1.*z27*z8*z85;
	z83 = -1.*z27*z8*z83;
	z176 = -1.*z27*z8*z81;
	z177 = z27*z8*z81;
	z178 = -1.*z27*z74*z8;
	z179 = -1.*z27*z51*z8;
	z8 = z27*z51*z8;
	z85 = -1.*z27*z85*z9;
	z81 = -1.*z27*z81*z9;
	z74 = -1.*z27*z74*z9;
	z180 = -1.*z27*z51*z9;
	z9 = z27*z51*z9;
	z27 = z106 + z50 + z52 + z56 + z77;
	z20 = z100 + z20 + z84 + z91;
	z50 = z107 + z52 + z59 + z78 + z92;
	z21 = z102 + z143 + z21 + z88;
	z51 = z104 + z144 + z55 + z90;
	z19 = z108 + z19 + z52 + z53 + z60;
	z10 = z10 + z13 + z146 + z23 + z96;
	z13 = z101 + z147 + z61 + z99;
	z53 = z103 + z109 + z67 + z86;
	z38 = z112 + z149 + z38;
	z55 = z113 + z151 + z64;
	z56 = z114 + z68 + z93;
	z30 = z101 + z117 + z154 + z30;
	z11 = z11 + z118 + z15 + z156 + z23;
	z15 = z105 + z119 + z66 + z89;
	z59 = z132 + z158 + z65;
	z41 = z133 + z160 + z41;
	z60 = z134 + z63 + z87;
	z61 = z103 + z138 + z163 + z70;
	z63 = z105 + z139 + z165 + z73;
	z12 = z12 + z140 + z17 + z23 + z82;
	z1 = z1 + z72 + z97;
	z2 = z2 + z58 + z98;
	z3 = z141 + z3 + z42;
	z17 = z106 + z128 + z145 + z24 + z44 + z48 + z49 + z52;
	z14 = z14 + z148 + z26 + z39 + z54 + z79;
	z23 = z115 + z128 + z49 + z57 + z80;
	z24 = z100 + z122 + z153 + z170 + z32 + z75;
	z31 = z123 + z157 + z171 + z31 + z45;
	z32 = z122 + z135 + z172 + z36;
	z25 = z102 + z124 + z162 + z173 + z25 + z76;
	z4 = z125 + z167 + z174 + z35 + z4;
	z7 = z124 + z142 + z43 + z7;
	z35 = z110 + z123 + z150 + z175 + z46;
	z36 = z107 + z116 + z129 + z155 + z49 + z52 + z62 + z83;
	z16 = z121 + z159 + z16 + z176 + z26 + z40;
	z39 = z120 + z129 + z177 + z49 + z69;
	z34 = z104 + z126 + z136 + z164 + z178 + z34;
	z5 = z127 + z168 + z179 + z37 + z5;
	z8 = z126 + z28 + z8 + z94;
	z28 = z111 + z125 + z152 + z47 + z85;
	z33 = z127 + z131 + z161 + z33 + z81;
	z22 = z108 + z130 + z137 + z166 + z22 + z49 + z52 + z74;
	z6 = z169 + z18 + z180 + z26 + z29 + z6;
	z9 = z130 + z49 + z71 + z9 + z95;

//	{{z27, z20, z21, z10, z30, z61, z38, z59, z1},
//	 {z20, z50, z51, z13, z11, z63, z55, z41, z2},
//	 {z21, z51, z19, z53, z15, z12, z56, z60, z3},
//   {z10, z13, z53, z17, z24, z25, z14, z31, z4},
//   {z30, z11, z15, z24, z36, z34, z35, z16, z5},
//   {z61, z63, z12, z25, z34, z22, z28, z33, z6},
//   {z38, z55, z56, z14, z35, z28, z23, z32, z7},
//   {z59, z41, z60, z31, z16, z33, z32, z39, z8},
//   {z1, z2, z3, z4, z5, z6, z7, z8, z9}}

	fK_3Body(0,0) = z27;
	fK_3Body(0,1) = z20;
	fK_3Body(0,2) = z21;
	fK_3Body(0,3) = z10;
	fK_3Body(0,4) = z30;
	fK_3Body(0,5) = z61;
	fK_3Body(0,6) = z38;
	fK_3Body(0,7) = z59;
	fK_3Body(0,8) = z1;
	fK_3Body(1,1) = z50;
	fK_3Body(1,2) = z51;
	fK_3Body(1,3) = z13;
	fK_3Body(1,4) = z11;
	fK_3Body(1,5) = z63;
	fK_3Body(1,6) = z55;
	fK_3Body(1,7) = z41;
	fK_3Body(1,8) = z2;
	fK_3Body(2,2) = z19;
	fK_3Body(2,3) = z53;
	fK_3Body(2,4) = z15;
	fK_3Body(2,5) = z12;
	fK_3Body(2,6) = z56;
	fK_3Body(2,7) = z60;
	fK_3Body(2,8) = z3;
	fK_3Body(3,3) = z17;
	fK_3Body(3,4) = z24;
	fK_3Body(3,5) = z25;
	fK_3Body(3,6) = z14;
	fK_3Body(3,7) = z31;
	fK_3Body(3,8) = z4;
	fK_3Body(4,4) = z36;
	fK_3Body(4,5) = z34;
	fK_3Body(4,6) = z35;
	fK_3Body(4,7) = z16;
	fK_3Body(4,8) = z5;
	fK_3Body(5,5) = z22;
	fK_3Body(5,6) = z28;
	fK_3Body(5,7) = z33;
	fK_3Body(5,8) = z6;
	fK_3Body(6,6) = z23;
	fK_3Body(6,7) = z32;
	fK_3Body(6,8) = z7;
	fK_3Body(7,7) = z39;
	fK_3Body(7,8) = z8;
	fK_3Body(8,8) = z9;
}

double SWDiamondT::Energy2Body(void)
{	
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;

	z1 = fLocd_2Body(0,0);
	z2 = fLocd_2Body(0,1);
	z3 = fLocd_2Body(0,2);
	z4 = fLocd_2Body(1,0);
	z5 = fLocd_2Body(1,1);
	z6 = fLocd_2Body(1,2);
	z7 = fLocX_2Body(0,0);
	z8 = fLocX_2Body(0,1);
	z9 = fLocX_2Body(0,2);
	z10 = fLocX_2Body(1,0);
	z11 = fLocX_2Body(1,1);
	z12 = fLocX_2Body(1,2);
	z10 = -z10;
	z11 = -z11;
	z12 = -z12;
	z4 = -z4;
	z5 = -z5;
	z6 = -z6;
	z1 = z1 + z10 + z4 + z7;
	z2 = z11 + z2 + z5 + z8;
	z3 = z12 + z3 + z6 + z9;
	z1 = z1*z1;
	z2 = z2*z2;
	z3 = z3*z3;
	z1 = z1 + z2 + z3;
	z1 = sqrt(z1);

	return( U2body(z1) );
}

void SWDiamondT::Force2Body(void)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;

	z1 = fLocd_2Body(0,0);
	z2 = fLocd_2Body(0,1);
	z3 = fLocd_2Body(0,2);
	z4 = fLocd_2Body(1,0);
	z5 = fLocd_2Body(1,1);
	z6 = fLocd_2Body(1,2);
	z7 = fLocX_2Body(0,0);
	z8 = fLocX_2Body(0,1);
	z9 = fLocX_2Body(0,2);
	z10 = fLocX_2Body(1,0);
	z11 = fLocX_2Body(1,1);
	z12 = fLocX_2Body(1,2);
	z10 = -z10;
	z11 = -z11;
	z12 = -z12;
	z4 = -z4;
	z5 = -z5;
	z6 = -z6;
	z1 = z1 + z10 + z4 + z7;
	z2 = z11 + z2 + z5 + z8;
	z3 = z12 + z3 + z6 + z9;
	z4 = z1*z1;
	z5 = z2*z2;
	z6 = z3*z3;
	z4 = z4 + z5 + z6;
	z4 = sqrt(z4);
	z5 = 1.0/z4;
	z4 = DU2body(z4);
	z6 = -z4*z5;
	z4 = z4*z5;
	z5 = z1*z6;
	z7 = z2*z6;
	z6 = z3*z6;
	z1 = z1*z4;
	z2 = z2*z4;
	z3 = z3*z4;
	
	//z1 = List(z5,z7,z6,z1,z2,z3);
	
	fF_2Body[0] = z5;
	fF_2Body[1] = z7;
	fF_2Body[2] = z6;
	fF_2Body[3] = z1;
	fF_2Body[4] = z2;
	fF_2Body[5] = z3;
}

void SWDiamondT::Stiffness2Body(void)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18, z19, z20, z21, z22, z23, z24;
	double z25, z26, z27;

	z1 = fLocd_2Body(0,0);
	z2 = fLocd_2Body(0,1);
	z3 = fLocd_2Body(0,2);
	z4 = fLocd_2Body(1,0);
	z5 = fLocd_2Body(1,1);
	z6 = fLocd_2Body(1,2);
	z7 = fLocX_2Body(0,0);
	z8 = fLocX_2Body(0,1);
	z9 = fLocX_2Body(0,2);
	z10 = fLocX_2Body(1,0);
	z11 = fLocX_2Body(1,1);
	z12 = fLocX_2Body(1,2);
	z10 = -z10;
	z11 = -z11;
	z12 = -z12;
	z4 = -z4;
	z5 = -z5;
	z6 = -z6;
	z1 = z1 + z10 + z4 + z7;
	z2 = z11 + z2 + z5 + z8;
	z3 = z12 + z3 + z6 + z9;
	z4 = z1*z1;
	z5 = z2*z2;
	z6 = z3*z3;
	z7 = z4 + z5 + z6;
	z8 = pow(z7,-1.5);
	z9 = pow(z7,-1.);
	z7 = sqrt(z7);
	z10 = 1.0/z7;
	z11 = DU2body(z7);
	z7  = DDU2body(z7);
	z8 = z11*z8;
	z12 = -z8;
	z10 = z10*z11;
	z13 = -z10;
	z14 = z1*z2*z7*z9;
	z11 = -z14;
	z16 = z1*z3*z7*z9;
	z15 = -z16;
	z18 = z2*z3*z7*z9;
	z17 = -z18;
	z20 = z4*z7*z9;
	z19 = -z20;
	z22 = z5*z7*z9;
	z21 = -z22;
	z7 = z6*z7*z9;
	z23 = -z7;
	z9 = z1*z12;
	z24 = z12*z2;
	z25 = z12*z4;
	z26 = z12*z5;
	z12 = z12*z6;
	z1 = z1*z8;
	z27 = z2*z8;
	z4 = z4*z8;
	z5 = z5*z8;
	z6 = z6*z8;
	z8 = z2*z9;
	z9 = z3*z9;
	z24 = z24*z3;
	z2 = z1*z2;
	z1 = z1*z3;
	z3 = z27*z3;
	z20 = z10 + z20 + z25;
	z22 = z10 + z22 + z26;
	z7 = z10 + z12 + z7;
	z4 = z13 + z19 + z4;
	z5 = z13 + z21 + z5;
	z6 = z13 + z23 + z6;
	z8 = z14 + z8;
	z9 = z16 + z9;
	z10 = z18 + z24;
	z2 = z11 + z2;
	z1 = z1 + z15;
	z3 = z17 + z3;
	
//	{{z20, z8, z9, z4, z2, z1},
//	 {z8, z22, z10, z2, z5, z3},
//	 {z9, z10, z7, z1, z3, z6},
//	 {z4, z2, z1, z20, z8, z9},
//	 {z2, z5, z3, z8, z22, z10},
//   {z1, z3, z6, z9, z10, z7}}
	
	fK_2Body(0,0) = z20;
	fK_2Body(0,1) = z8;
	fK_2Body(0,2) = z9;
	fK_2Body(0,3) = z4;
	fK_2Body(0,4) = z2;
	fK_2Body(0,5) = z1;
	fK_2Body(1,1) = z22;
	fK_2Body(1,2) = z10;
	fK_2Body(1,3) = z2;
	fK_2Body(1,4) = z5;
	fK_2Body(1,5) = z3;
	fK_2Body(2,2) = z7;
	fK_2Body(2,3) = z1;
	fK_2Body(2,4) = z3;
	fK_2Body(2,5) = z6;
	fK_2Body(3,3) = z20;
	fK_2Body(3,4) = z8;
	fK_2Body(3,5) = z9;
	fK_2Body(4,4) = z22;
	fK_2Body(4,5) = z10;
	fK_2Body(5,5) = z7;
}

/* print connectivity element data */
void SWDiamondT::PrintConnectivityData(ostream& out)
{
	out << "\n Number of 3 body interactions . . . . . . . . . = " << fNodes_3Body.MajorDim() << '\n';
	out <<   " Number of 2 body interactions . . . . . . . . . = " << fNodes_2Body.MajorDim() << '\n';

	/* 3-body connectivities */
	out << "\n 3-body Connectivities:\n\n";
	out << setw(kIntWidth) << "no.";
	for (int i = 1; i <= 3; i++)
	{
		out << setw(kIntWidth - 2) << "n[";
		out << i << "]";
	}
	out << '\n';
	fNodes_3Body.WriteNumbered(out);

	/* 2-body connectivities */
	out << "\n 2-body Connectivities:\n\n";
	out << setw(kIntWidth) << "no.";
	for (int j = 1; j <= 2; j++)
	{
		out << setw(kIntWidth - 2) << "n[";
		out << j << "]";
	}
	out << '\n';
	fNodes_2Body.WriteNumbered(out);
}	

/***********************************************************************
* Private
***********************************************************************/

/* 3 body potentials */
double SWDiamondT::U3body(double r1, double r2, double c12) const
{
	return( pow(1.0/3.0 + c12, 2.0)*feps*flambda*
	         ( (r1 < frcut*fa) ? exp(fgamma/(-frcut + r1/fa)) : 0.0 )*
( (r2 < frcut*fa) ? exp(fgamma/(-frcut + r2/fa)) : 0.0) );
}
	
	/* 1st derivs */
double SWDiamondT::Dr1U3body(double r1, double r2, double c12) const
{
	return( pow(1.0 + 3.0*c12, 2.0)*feps*flambda*
	        ( (r1 < frcut*fa) ? -(fa*fgamma*exp(fgamma/(-frcut + r1/fa))/
	                          pow(-(fa*frcut) + r1, 2.0)) : 0.0)*
( (r2 < frcut*fa) ? exp(fgamma/(-frcut + r2/fa)) : 0.0 )/9.0 );
}

double SWDiamondT::Dr2U3body(double r1, double r2, double c12) const
{
	return( pow(1.0 + 3.0*c12, 2.0)*feps*flambda*
	         ( (r1 < frcut*fa) ? exp(fgamma/(-frcut + r1/fa)) : 0.0 )*
( (r2 < frcut*fa) ? -(fa*fgamma*exp(fgamma/(-frcut + r2/fa))/
pow(-(fa*frcut) + r2, 2.0)) : 0.0)/9.0 );
}

double SWDiamondT::Dc12U3body(double r1, double r2, double c12) const
{
	return( 2.0*(1.0/3.0 + c12)*feps*flambda*
	        ( (r1 < frcut*fa) ? exp(fgamma/(-frcut + r1/fa)) : 0.0 )*
( (r2 < frcut*fa) ? exp(fgamma/(-frcut + r2/fa)) : 0.0 ) );
}

	/* 2nd derivs */
double SWDiamondT::DDr1U3body(double r1, double r2, double c12) const
{
	return( pow(1.0 + 3.0*c12, 2.0)*feps*flambda*
	         ( (r1 < frcut*fa) ? fa*fgamma*(fa*fgamma - 2*fa*frcut + 2*r1)*
	            exp(fgamma/(-frcut + r1/fa))/pow(-(fa*frcut) + r1, 4.0) : 0.0)*
	         ( (r2 < frcut*fa) ? exp(fgamma/(-frcut + r2/fa)) : 0.0)/9.0 );
}

double SWDiamondT::DDr2U3body(double r1, double r2, double c12) const
{
	return( pow(1.0 + 3.0*c12, 2.0)*feps*flambda*
	         ( (r1 < frcut*fa) ? exp(fgamma/(-frcut + r1/fa)) : 0.0 )*
( (r2 < frcut*fa) ? fa*fgamma*(fa*fgamma - 2*fa*frcut + 2*r2)*
exp(fgamma/(-frcut + r2/fa))/pow(-(fa*frcut) + r2, 4.0) : 0.0 )/9.0 );
}

double SWDiamondT::DDc12U3body(double r1, double r2, double c12) const
{
#pragma unused(c12) //cos_12 doesn't appear in second derivative

	return( 2.0*feps*flambda*
	        ( (r1 < frcut*fa) ? exp(fgamma/(-frcut + r1/fa)) : 0.0 )*
( (r2 < frcut*fa) ? exp(fgamma/(-frcut + r2/fa)) : 0.0 ) );
}
		
	/* mixed derivs */
double SWDiamondT::Dr1Dr2U3body(double r1, double r2, double c12) const
{
	return( pow(1.0 + 3.0*c12, 2.0)*feps*flambda*
	         ( (r1 < frcut*fa) ? -(fa*fgamma*exp(fgamma/(-frcut + r1/fa))/
	                            pow(-(fa*frcut) + r1,2)) : 0.0 )*
( (r2 < frcut*fa) ? -(fa*fgamma*exp(fgamma/(-frcut + r2/fa))/
pow(-(fa*frcut) + r2,2)) : 0.0 )/9.0 );
}

double SWDiamondT::Dr1Dc12U3body(double r1, double r2, double c12) const
{
	return( 2.0*(1.0/3.0 + c12)*feps*flambda*
	         ( (r1 < frcut*fa) ? -(fa*fgamma*exp(fgamma/(-frcut + r1/fa))/
	                            pow(-(fa*frcut) + r1, 2.0)) : 0.0 )*
( (r2 < frcut*fa) ? exp(fgamma/(-frcut + r2/fa)) : 0.0 ) );
}

double SWDiamondT::Dr2Dc12U3body(double r1, double r2, double c12) const
{
	return( 2.0*(1.0/3.0 + c12)*feps*flambda*
	        ( (r1 < frcut*fa) ? exp(fgamma/(-frcut + r1/fa)) : 0.0 )*
( (r2 < frcut*fa) ? -(fa*fgamma*exp(fgamma/(-frcut + r2/fa))/
pow(-(fa*frcut) + r2, 2.0)) : 0.0 ) );
}

/* compute entire Hessian at once */
void SWDiamondT::ComputeHessian_3Body(dMatrixT& hessian,
	double r1, double r2, double c12)
{
	double z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12;
	double z13, z14, z15, z16, z17, z18;

	z1 = 0.3333333333333333 + c12;
	z2 = pow(fa,-2.);
	z3 = pow(fa,-1.);
	z4 = pow(fgamma,2.);
	z5 = -frcut;
	z6 = pow(z1,2.);
	z7 = r1*z3;
	z8 = r2*z3;
	z7 = z5 + z7;
	z5 = z5 + z8;
	z8 = pow(z5,-4.);
	z9 = pow(z5,-3.);
	z10 = pow(z5,-2.);
	z5 = pow(z5,-1.);
	z11 = pow(z7,-4.);
	z12 = pow(z7,-3.);
	z13 = pow(z7,-2.);
	z7 = pow(z7,-1.);
	z5 = fgamma*z5;
	z7 = fgamma*z7;
	z5 = exp(z5);
	z7 = exp(z7);
	z14 = (r1 < frcut) ? z7 : 0.0;
	z15 = (r2 < frcut) ? z5 : 0.0;
	z16 = 2.*fgamma*z2;
	z17 = z2*z5;
	z18 = 2.*feps*flambda*z14*z15;
	z9 = z16*z5*z9;
	z12 = z12*z16*z7;
	z2 = z11*z2*z4*z7;
	z7 = (r1 < frcut) ? -fgamma*z13*z3*z7 : 0.0;
	z3 = (r2 < frcut) ? -fgamma*z10*z3*z5 : 0.0;
	z1 = 2.*feps*flambda*z1;
	z5 = z1*z15*z7;
	z1 = z1*z14*z3;
	z4 = z17*z4*z8;
	z3 = feps*flambda*z3*z6*z7;
	z2 = (r1 < frcut) ? z12 + z2 : 0;
	z2 = feps*flambda*z15*z2*z6;
	z4 = (r2 < frcut) ? z4 + z9 : 0;

	//result:
	//{{z2, z3, z5},
	// {z3, feps flambda z14 z4 z6, z1},
	// {z5, z1, z18}}

	/* assign upper triangle only */
	hessian(0,0) = z2;
	hessian(1,1) = feps*flambda*z14*z4*z6;
	hessian(2,2) = z18;
	hessian(1,2) = z1;
	hessian(0,2) = z5;
	hessian(0,1) = z3;

}	

	/* 2 body potential and derivatives */
double SWDiamondT::U2body(double r) const
{
	return( fA*feps*(-1.0 + pow(fa,4.0)*fB/pow(r,4.0))*
	         ( (r < frcut*fa) ? exp(fdelta/(-frcut + r/fa)) : 0.0 ) );
}

double SWDiamondT::DU2body(double r) const
{
	return( -4.0*pow(fa,4.0)*fA*fB*feps*
	          ( (r < frcut*fa) ? exp(fdelta/(-frcut + r/fa)) : 0.0 )/pow(r,5.0) +
fA*feps*(-1.0 + pow(fa,4.0)*fB/pow(r,4.0))*
( (r < frcut*fa) ? -(fa*fdelta*exp(fdelta/(-frcut + r/fa))/
pow(-(fa*frcut) + r,2)) : 0.0) );
}

double SWDiamondT::DDU2body(double r) const
{
	return( 20.0*pow(fa,4.0)*fA*fB*feps*
	         ( (r < frcut*fa) ? exp(fdelta/(-frcut + r/fa)) : 0.0)/pow(r,6.0) -
8.0*pow(fa,4.0)*fA*fB*feps*
( (r < frcut*fa) ? -(fa*fdelta*exp(fdelta/(-frcut + r/fa))/
pow(-(fa*frcut) + r,2.0)) : 0.0 )/pow(r,5.0) +
fA*feps*(-1.0 + pow(fa,4.0)*fB/pow(r,4.0))*
( (r < frcut*fa) ? fa*fdelta*(fa*fdelta - 2.0*fa*frcut + 2.0*r)*
exp(fdelta/(-frcut + r/fa))/
pow(-(fa*frcut) + r,4.0): 0.0 ) );
}
