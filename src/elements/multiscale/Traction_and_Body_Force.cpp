/* $Id: Traction_and_Body_Force.cpp,v 1.4 2004-07-15 08:28:27 paklein Exp $ */
#include "StaggeredMultiScaleT.h"

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "ShapeFunctionT.h"
#include "Traction_CardT.h"
#include "iAutoArrayT.h"
#include "ScheduleT.h"

#include "VariArrayT.h"
#include "nVariArray2DT.h"
#include "VariLocalArrayT.h"

//##################################################################################
//###### Traciton B.C. Methods (Cut and Paste from StaggeredMultiScaleT) ##############
//##################################################################################

//---------------------------------------------------------------------

void StaggeredMultiScaleT::EchoTractionBC(ifstreamT& in, ostream& out)
{
ExceptionT::GeneralFail("StaggeredMultiScaleT::EchoTractionBC", "out of date");
#if 0	
	const char caller[] = "StaggeredMultiScaleT::EchoTractionBC";
	
	/* read data from parameter file */
	int numlines, numsets;
	ModelManagerT& model = ElementSupport().Model();
	model.ReadNumTractionLines (in, numlines, numsets);

	if (numlines > 0)
	  {
	    /* temp space */
	    ArrayT<StringT> block_ID(numlines);
	    ArrayT<iArray2DT> localsides (numlines);
	    iArrayT LTf (numlines);
	    ArrayT<Traction_CardT::CoordSystemT> coord_sys (numlines);
	    ArrayT<dArray2DT> values (numlines);

	    /* nodes on element facets */
	    iArrayT num_facet_nodes;
	    fShapes->NumNodesOnFacets (num_facet_nodes);
	    
	    /* read data by blocks */
	    int line = 0;
	    int count = 0;
	    for (int blockset = 0; blockset < numsets; blockset++)
	      {
		/* read num of cards in each block */
		int setsize = -1;
		StringT set_ID;
		model.ReadTractionSetData (in, set_ID, setsize);
		for (int card=0; card < setsize; card++)
		  {
		    /* read side set for that card */
		    block_ID[line] = set_ID;
		    model.ReadTractionSideSet (in, block_ID[line], localsides[line]);

		    /* increment count */
		    int num_sides = localsides[line].MajorDim();
		    count += num_sides;

		    /* read data for that card */
		    in >> LTf[line] >> coord_sys[line];

		    /* skip if empty */
		    int num_nodes;
		    if (num_sides > 0)
		      {
			iArray2DT& side_set = localsides[line];

			/* switch to group numbering */
			const ElementBlockDataT& block_data = BlockData (block_ID[line]);
			iArrayT elems (num_sides);
			side_set.ColumnCopy (0, elems);

			/* check */
			int min, max;
			elems.MinMax (min, max);
			if (min < 0 || max > block_data.Dimension())
				ExceptionT::BadInputValue(caller, 
					"node numbers {%d,%d} are out of range in dataline %d",
					min, max, line);

			/* shift */
			elems += block_data.StartNumber();
			side_set.SetColumn (0, elems);

			/* all facets in set must have the same number of nodes */
			num_nodes = num_facet_nodes [side_set (0,1)];
			for (int f=0; f < num_sides; f++)
			  if (num_facet_nodes[side_set(f,1)] != num_nodes)
			  	ExceptionT::BadInputValue(caller,
			  		"sides specified in line %d have differing numbers of nodes", line);
		      }
		    else
		      {
			/* still check number of facet nodes */
			int min, max;
			num_facet_nodes.MinMax (min, max);
			if (min != max)
				ExceptionT::BadInputValue(caller,
					"cannot determine number of facet nodes for empty side set at line %d", line);
			else
			  num_nodes = min;
		      }

		    /* read traction values */
		    dArray2DT& valueT = values[line];
		    valueT.Allocate (num_nodes, NumDOF());
		    in >> valueT;
		    //NOTE - cannot simply clear to the end of the line with empty side sets
		    //       because the tractions may be on multiple lines

		    line++;
		  }
	      }

	    /* allocate all traction BC cards */
	    fTractionList.Allocate (count);

	    /* correct numbering offset */
	    LTf--;

	    if (count > 0)
	      {
		out << '\n';
		out << setw (kIntWidth) << "no.";
		fTractionList[0].WriteHeader (out, NumDOF());

		iArrayT loc_node_nums;
		int dex = 0;
		for (int ii=0; ii < numlines; ii++)
		  {
		    /* set traction BC cards */
		    iArray2DT& side_set = localsides[ii];
		    int numsides = side_set.MajorDim();
		    for (int j=0; j < numsides; j++)
		      {
			out << setw (kIntWidth) << j+1;

			/* get facet local node numbers */
			fShapes->NodesOnFacet (side_set (j, 1), loc_node_nums);

			/* set and echo */
			fTractionList[dex++].EchoValues (ElementSupport(), side_set(j,0), side_set (j,1), LTf[ii],
							 coord_sys[ii], loc_node_nums, values[ii], out);
		      }
		    out << endl;
		  }
	      }
	  }


	if (NumSD() != NumDOF())
	{
		/* check coordinate system specifications */
		for (int i = 0; i < fTractionList.Length(); i++)
			if (fTractionList[i].CoordSystem() != Traction_CardT::kCartesian)
			{
				cout << "\n StaggeredMultiScaleT::EchoTractionBC: coordinate system must be\n"
				     <<   "    Cartesian:" << Traction_CardT::kCartesian
				     << " if (spatial dimensions != degrees of freedom)\n"
				     <<   "    for card " << i+1 << endl;
				throw ExceptionT::kBadInputValue;
			}
	}
#endif
}

//---------------------------------------------------------------------

/* update traction BC data */
void StaggeredMultiScaleT::SetTractionBC(void)
{
//NOTE: With the possibility of variable global node numbers and
//		and equations, we assume as little as possible here with
//      regard to the validity of the node/equation numbers, requiring
//      only that NodesX in the element cards has the correct global
//      node numbers.

	/* dimensions */
	int ndof = NumDOF();

	/* echo values */
	iArray2DT nd_tmp, eq_tmp;
	for (int i = 0; i < fTractionList.Length(); i++)
	{
		Traction_CardT& BC_card = fTractionList[i];
			
		/* traction element/facet */
		int elem, facet;
		BC_card.Destination(elem, facet);

		/* set global node numbers */
		const iArrayT& loc_nodes = BC_card.LocalNodeNumbers();
		int nnd = loc_nodes.Length();
		
		iArrayT& nodes = BC_card.Nodes();
		nodes.Dimension(nnd);
		nodes.Collect(loc_nodes, fElementCards[elem].NodesX());
		
		/* set global equation numbers */
		iArrayT& eqnos = BC_card.Eqnos();
		eqnos.Dimension(ndof*nnd);
		
		/* get from node manager */
		nd_tmp.Set(1, nnd, nodes.Pointer());
		eq_tmp.Set(1, ndof*nnd, eqnos.Pointer());
		fCoarse.SetLocalEqnos(nd_tmp, eq_tmp);
	}

	/* set flag */
	fTractionBCSet = 1;
}

//---------------------------------------------------------------------

/* compute contribution to RHS from traction BC's */
void StaggeredMultiScaleT::ApplyTractionBC(void)
{
	if (fTractionList.Length() > 0)
	{
		/* dimensions */
		int nsd = NumSD();
		int ndof = NumDOF();
	
		/* update equation numbers */
		if (!fTractionBCSet) SetTractionBC();
	
		/* force vector */
		dArrayT rhs;
		VariArrayT<double> rhs_man(25, rhs);
		
		/* local coordinates */
		LocalArrayT coords(LocalArrayT::kInitCoords);
		VariLocalArrayT coord_man(25, coords, nsd);
		ElementSupport().RegisterCoordinates(coords);
		
		/* nodal tractions */
		LocalArrayT tract(LocalArrayT::kUnspecified);
		VariLocalArrayT tract_man(25, tract, ndof);

		/* integration point tractions */
		dArray2DT ip_tract;
		nVariArray2DT<double> ip_tract_man(25, ip_tract, ndof);
		dArrayT tract_loc, tract_glb(ndof);
		dMatrixT Q(ndof);
		
		/* Jacobian of the surface mapping */
		dMatrixT jacobian(nsd, nsd-1);
		
		for (int i = 0; i < fTractionList.Length(); i++)
		{
			const Traction_CardT& BC_card = fTractionList[i];

			/* dimension */
			const iArrayT& nodes = BC_card.Nodes();
			int nnd = nodes.Length();
			rhs_man.SetLength(nnd*ndof, false);
			coord_man.SetNumberOfNodes(nnd);
			tract_man.SetNumberOfNodes(nnd);
			
			/* local coordinates */
			coords.SetLocal(nodes);

			/* nodal traction vectors: (ndof x nnd) */
			BC_card.CurrentValue(tract);
			
			/* BC destination */
			int elem, facet;
			BC_card.Destination(elem, facet);
			
			/* default thickness */
			double thick = 1.0;
			
			/* boundary shape functions */
			const ParentDomainT& surf_shape = ShapeFunction().FacetShapeFunction(facet);
			int nip = surf_shape.NumIP();
			
			/* all ip tractions: (nip x ndof) */
			ip_tract_man.SetMajorDimension(nip, false);
			surf_shape.Interpolate(tract, ip_tract);

			/* traction vector coordinate system */
			if (BC_card.CoordSystem() == Traction_CardT::kCartesian)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian);
	
					/* ip weight */
					double jwt = detj*w[j]*thick;
					
					/* ip traction */
					const double* tj = ip_tract(j);
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else if (BC_card.CoordSystem() == Traction_CardT::kLocal)
			{
				/* integrate */			
				rhs = 0.0;
				const double* w = surf_shape.Weight();
				for (int j = 0; j < nip; j++)
				{
					/* coordinate mapping */
					surf_shape.DomainJacobian(coords, j, jacobian);
					double detj = surf_shape.SurfaceJacobian(jacobian, Q);
	
					/* ip weight */
					double jwt = detj*w[j]*thick;
					
					/* transform ip traction out of local frame */
					ip_tract.RowAlias(j, tract_loc);
					Q.Multx(tract_loc, tract_glb);

					/* ip traction */
					const double* tj = tract_glb.Pointer();
					
					/* accumulate */
					for (int l = 0; l < ndof; l++)
					{
						/* nodal shape function */
						const double* Na = surf_shape.Shape(j);
					
						double* prhs = rhs.Pointer(l);
						double  fact = jwt*(*tj++);
						for (int k = 0; k < nnd; k++)
						{
							*prhs += fact*(*Na++);
							prhs += ndof;
						}
					}				
				}
			}
			else
				throw ExceptionT::kGeneralFail;

			/* assemble into coarse scale equations */
			ElementSupport().AssembleRHS(fCoarse.Group(), rhs, BC_card.Eqnos());
		}
	}
}


//##################################################################################
//###### Body Force Methods (Cut and Paste from StaggeredMultiScaleT) #################
//##################################################################################
//---------------------------------------------------------------------

void StaggeredMultiScaleT::EchoBodyForce(ifstreamT& in, ostream& out)
{
	const char caller[] = "StaggeredMultiScaleT::EchoBodyForce";

	/* schedule number and body force vector */
	int n_sched;
	in >> n_sched >> fBody;		
	n_sched--;

	/* no LTf => no body force */
	if (n_sched < 0) 
		fBody = 0.0;
	else
	{
		fBodySchedule = ElementSupport().Schedule(n_sched);
		if (!fBodySchedule)
			ExceptionT::BadInputValue(caller, "could not resolve schedule %d", n_sched + 1);
	}
	
	out << "\n Body force vector:\n";
	out << " Body force load-time function number. . . . . . = " << n_sched + 1<< '\n';
	out << " Body force vector components:\n";
	for (int j = 0 ; j < NumDOF(); j++)
	{
		out << "   x[" << j+1 << "] direction. . . . . . . . . . . . . . . . = ";
		out << fBody[j] << '\n';
	}
	out.flush();   	   	
}

//---------------------------------------------------------------------
// Purpose:	body force data is extracted from global data storage and put in 
// 					LocalArrayT for current element 
// Input: 	none 
// Output: 	body_force (n_en x n_dof) (i.e at nodes)
// Action:	global data from fBodySchedule put into LocalArrayT

/* add contribution from the body force */
void StaggeredMultiScaleT::AddBodyForce(LocalArrayT& body_force) const
{
	if (fBodySchedule)
	{
		int ndof = NumDOF();
		int nen = body_force.NumberOfNodes();
		double loadfactor = fBodySchedule->Value();
		double* p = body_force.Pointer();

		for (int i = 0; i < ndof; i++)
		{
			double temp = -fBody[i]*loadfactor;
			for (int j = 0; j < nen; j++)
				*p++ = temp;
		}
	}
}

//---------------------------------------------------------------------
//
// Input: 	body forces either at nodes or at ip.
// Output:  none
// Action:  fRHS is updated to include body forces

/* calculate the body force contribution (constM must be the density rho) */
void StaggeredMultiScaleT::FormMa(MassTypeT mass_type, double constM, 
	const LocalArrayT* nodal_values,
	const dArray2DT* ip_values)
{
	const char caller[] = "StaggeredMultiScaleT::FormMa";

	/* quick exit */
	if (!nodal_values && !ip_values) return;

#if __option(extended_errorcheck)
	/* dimension checks */
	if (nodal_values && 
		fRHS.Length() != nodal_values->Length()) 
		ExceptionT::SizeMismatch(caller);

	if (ip_values &&
		(ip_values->MajorDim() != fShapes->NumIP() ||
		 ip_values->MinorDim() != NumDOF()))
		ExceptionT::SizeMismatch(caller);
#endif

	switch (mass_type)
	{
		case kConsistentMass:
		{
			int ndof = NumDOF();
			int  nen = NumElementNodes();

			const double* Det    = fShapes->IPDets();
			const double* Weight = fShapes->IPWeights();

			fShapes->TopIP();
			while (fShapes->NextIP())
			{					
				/* interpolate nodal values to ip */
				if (nodal_values)
					fShapes->InterpolateU(*nodal_values, fDOFvec);
					
				/* ip sources */
				if (ip_values)
					fDOFvec -= (*ip_values)(fShapes->CurrIP());

				/* accumulate in element residual force vector */				
				double*	res      = fRHS.Pointer();
				const double* Na = fShapes->IPShapeU();
				
				double temp = constM*(*Weight++)*(*Det++);				
				for (int lnd = 0; lnd < nen; lnd++)
				{
					double temp2 = temp*(*Na++);
					double* pacc = fDOFvec.Pointer();

					for (int dof = 0; dof < ndof; dof++)			
						*res++ += temp2*(*pacc++);
				}
			}
			break;
		}	
		case kLumpedMass:
		{
			ExceptionT::GeneralFail(caller, "no lumped mass");
			break;
		}
	}
}


//------ End 
