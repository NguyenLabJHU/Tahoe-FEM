/* $Id: FEManagerT_THK.cpp,v 1.7 2003-07-11 16:45:19 hspark Exp $ */
#include "FEManagerT_THK.h"
#ifdef BRIDGING_ELEMENT

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "FieldT.h"
#include "StringT.h"
#include "ParticlePairT.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "RaggedArray2DT.h"

using namespace Tahoe;

const double tol = 1.0e-8;   // for neighbor searching tolerance
const double root32 = sqrt(3.0)/2.0;    // for neighbor searching tolerance

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, bridging_input),
	fTheta(2)
{

}

/* initialize members */
void FEManagerT_THK::Initialize(InitCodeT init)
{
	/* inherited */
	FEManagerT_bridging::Initialize(init);
	
	// read other parameters and initialize data
	ifstreamT& in = Input();
	in >> fNcrit;
	int num_neighbors = 2 * fNcrit + 1;   // maximum number of neighbors per atom in 2D
		  
	/* obtain list of atoms on which BC's will be applied in FEManagerT_THK */
	ArrayT<StringT> id_list;
        
	ModelManagerT* model = FEManagerT::ModelManager();
        
	/* read node set indexes */
	model->NodeSetList(in, id_list);
	int numsets = id_list.Length();		// number of MD THK boundary node sets
 
	/* collect sets - currently assuming exactly 2 MD THK boundaries in 2D */
	/* further assumes bottom is first node set, top is second node set */
	fBottomatoms = model->NodeSet(id_list[0]);
	fTopatoms = model->NodeSet(id_list[1]);
	fBottomrow.Dimension(fBottomatoms.Length());
	fBottomrow.SetValueToPosition();
	fToprow.Dimension(fTopatoms.Length());
	fToprow.SetValueToPosition();
	fToprow+=fBottomrow.Length();
	int num_nodes_top = fTopatoms.Length();		// total number of boundary atoms on top layer
	int num_nodes_bottom = fBottomatoms.Length();  // total number of boundary atoms on bottom layer
	fBottom.Dimension(num_nodes_bottom, num_neighbors);
	fBottom = -1;  // -1 -> no neighbor
	fTop.Dimension(num_nodes_top, num_neighbors);
	fTop = -1;  // -1 -> no neighbor

	double lparam, xbottom, xtop, ybottom, ytop;
	in >> lparam;    // Read lattice parameter in from input file
	
	NodeManagerT* node = FEManagerT::NodeManager();
	const dArray2DT& initcoords = node->InitialCoordinates();  // init coords for all atoms
	dArrayT bacoordst, bacoordsb, currcoordst, currcoordsb;   // boundary atom coordinates
	
	/* now loop over boundary atoms to find neighbors - specialized for 2D FCC lattices */
	/* assume top and bottom have same number of atoms */
	for (int i = 0; i < fBottomatoms.Length(); i++)	// changed from fBottom.MajorDim()
	{
		/* obtain coordinates of each bottom and top boundary atom */
		initcoords.RowAlias(fBottomatoms[i], bacoordsb);
		initcoords.RowAlias(fTopatoms[i], bacoordst);
		
		for (int j = -fNcrit; j <= fNcrit; j++)
		{
			if ((i+j >= 0) && (i+j < fBottomatoms.Length()))	// check with boundary only
			{	
				initcoords.RowAlias(fBottomatoms[i+j], currcoordsb);
				initcoords.RowAlias(fTopatoms[i+j], currcoordst);
				fBottom(i,j+fNcrit) = fBottomatoms[i+j];	// +1
				fTop(i,j+fNcrit) = fTopatoms[i+j];	// +1
			}
		
		
		}
	}
	
	/* compute theta tables */
	StringT path;
	path.FilePath(in.filename());
	StringT data_file;
	in >> data_file;
	data_file.ToNativePathName();
	data_file.Prepend(path);
	ComputeThetaTables(data_file);
}

/* initialize the current time increment for all groups */
ExceptionT::CodeT FEManagerT_THK::InitStep(void)
{
	/* inherited */
	ExceptionT::CodeT result = FEManagerT_bridging::InitStep();
	return result;
}

/* close the current time increment for all groups */
ExceptionT::CodeT FEManagerT_THK::CloseStep(void)
{
	/* inherited */
	ExceptionT::CodeT result = FEManagerT_bridging::CloseStep();
	return result;
}

/* return iArrayT of boundary and ghost atom numbers */
const iArrayT& FEManagerT_THK::InterpolationNodes(void)
{
	const iArrayT& ghost = GhostNodes();
	
	fInterpolationNodes.Dimension(ghost.Length()+fBottomatoms.Length()+fTopatoms.Length());
	fInterpolationNodes.CopyIn(0,ghost);	// copy in ghost atoms first
	fInterpolationNodes.CopyIn(ghost.Length(), fBottomatoms);	// then copy in bottom atoms
	fInterpolationNodes.CopyIn(ghost.Length()+fBottomatoms.Length(), fTopatoms);	// finally top atoms
	
	return fInterpolationNodes;
}

/* predictor and corrector routine for FEM solution interpolated to MD boundary atoms.  both predictor
   and corrector are done in one step due to constant coarse scale acceleration assumption.  */
void FEManagerT_THK::BAPredictAndCorrect(double timestep, dArray2DT& badisp, dArray2DT& bavel, dArray2DT& baacc)
{
	/* displacement predictor (and corrector) */
	badisp.AddCombination(timestep, bavel, .5*timestep*timestep, baacc);
	
	/* velocity predictor (and corrector) */
	bavel.AddScaled(timestep, baacc);
}

/* calculate external force on MD boundary atoms using time history kernel approach */
const dArray2DT& FEManagerT_THK::THKForce(const dArray2DT& badisp)
{
	StringT bridging_field = "displacement";
	
	/* badisp is in format bottom displacements first, then top */
	fTHKforce.Dimension(badisp.MajorDim(), 2);  

	/* sort top and bottom displacements into separate arrays - may not be necessary */
	dArray2DT topdisp, bottomdisp;
	dArrayT atomdisp, femdisp(2), diff(2);
	topdisp.Dimension(fTopatoms.Length(), 2);
	bottomdisp.Dimension(fBottomatoms.Length(), 2);
	bottomdisp.RowCollect(fBottomrow, badisp);
	topdisp.RowCollect(fToprow, badisp);

	/* access the actual MD displacements */
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);
	dArray2DT mddisp = (*atomfield)[0];

	const int stepnum = FEManagerT::StepNumber();  // to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();  // delta t_md

	/* loop over boundary atoms first, store q - ubar */
	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, 2);	// copy all rows of history except last row

	/* Calculate q - ubar */
	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		if (stepnum < fNumstep_crit)   
		{
			/* bottom layer of atoms */
			mddisp.RowAlias(fBottomatoms[i], atomdisp);
			bottomdisp.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTableb[i].SetRow(stepnum, diff);
	
			/* top layer of atoms */
			mddisp.RowAlias(fTopatoms[i], atomdisp);   
			topdisp.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTablet[i].SetRow(stepnum, diff);
		}
		else	// time greater than critical normalized time
		{
			/* bottom layer of atoms */
			mddisp.RowAlias(fBottomatoms[i], atomdisp);
			bottomdisp.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTableb[i]);
			fHistoryTableb[i].BlockRowCopyAt(shift, 0);
			fHistoryTableb[i].SetRow(fNumstep_crit-1, diff);
	
			/* top layer of atoms */
			mddisp.RowAlias(fTopatoms[i], atomdisp);
			topdisp.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTablet[i]);
			fHistoryTablet[i].BlockRowCopyAt(shift, 0);
			fHistoryTablet[i].SetRow(fNumstep_crit-1, diff);
		}
	}

	dMatrixT theta(2);
	dArrayT force1(2), force2(2), force3b(2), force3t(2);

	/* loop over boundary atoms, calculate THK force, i.e. Theta(t-tau)*(q-ubar)(tau)*(delta tau) */
	for (int i = 0; i < fBottom.MajorDim(); i++)
	{
		force3b = 0.0;
		force3t = 0.0;
		for (int j = -fNcrit; j <= fNcrit; j++)
		{
			if (fBottom(i,j+fNcrit) == -1)
				int blah = 0;	// do nothing
			else
			{
				const dArray2DT& theta_b = fThetaTable[j+fNcrit];
				const dArray2DT& disp_b = fHistoryTableb[i+j];
				const dArray2DT& disp_t = fHistoryTablet[i+j];
				const dArray2DT& theta_t = fThetaTable[fNcrit-j];  
		
				/* calculate THK force using theta here */
				if (stepnum < fNumstep_crit) 
				{
					for (int k = 0; k < stepnum; k++)
					{
						/* bottom row */
						theta.Set(2,2,theta_b(k));
						disp_b.RowAlias(stepnum-k, force1);
						theta.Multx(force1, force2);
						force3b.AddScaled(timestep, force2);
					
						/* top row */
						theta.Set(2,2,theta_t(k));
						disp_t.RowAlias(stepnum-k, force1);  
						theta.Multx(force1, force2);  
						force3t.AddScaled(timestep, force2);   
					}
				}
				else	// normalized time greater than critical value
				{
					for (int k = 0; k < fNumstep_crit; k++)
					{
						/* bottom row */
						theta.Set(2,2,theta_b(k));
						disp_b.RowAlias(fNumstep_crit-k-1, force1);
						theta.Multx(force1, force2);
						force3b.AddScaled(timestep, force2);
			
						/* top row */
						theta.Set(2,2,theta_t(k));
						disp_t.RowAlias(fNumstep_crit-k-1, force1);
						theta.Multx(force1, force2);
						force3t.AddScaled(timestep, force2);
					}
				}
			}
		}
		/* add THK force to fTHKForce - bottom atoms */
		fTHKforce.SetRow(i, force3b);
		
		/* add THK force to fTHKForce - top atoms */
		fTHKforce.SetRow(i+fBottomatoms.Length(), force3t);
	}
	return fTHKforce;
}

/* function to modify initial positions via Gaussian displacement in y-direction */
const dArray2DT& FEManagerT_THK::GaussianWave(void) 
{
        double sig = 15.0;
        double H = sig/4.0;
        double A = 1.5e-1;
        double b = 0.1;
        double Lc = 5.0 * sig;
        double uc = A * exp(-pow(Lc/sig,2));
        double pi = acos(-1.0);
        double xi, disp;
	
        /* loop over all atoms, apply initial displacement */
        /* initial displacement identical to Wagner and Liu, but applied in y-direction only on all atoms */
        NodeManagerT* node = FEManagerT::NodeManager();
        const dArray2DT& initcoords = node->InitialCoordinates();  // init coords for all atoms
        fGaussdisp.Dimension(initcoords.MajorDim(), initcoords.MinorDim());
        fGaussdisp = 0.0;
	
	for (int i = 0; i < initcoords.MajorDim(); i++)
	{
		xi = fabs(initcoords(i,1));	// check y-coords only
		if (xi <= Lc)
		{
			disp = A/(A-uc)*(1.0+b*cos(2*pi*xi/H))*(A*exp(-pow(xi/sig,2))-uc);
			fGaussdisp(i,1) = disp;
		}
	}
	return fGaussdisp;
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute theta tables */
void FEManagerT_THK::ComputeThetaTables(const StringT& data_file)
{
	const char caller[] = "FEManagerT_THK::ComputeThetaTables";
	ifstreamT data(data_file);
	if (!data.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
		
	/* dimensions */
	double pi = acos(-1.0);
	int n_neighbor = fNcrit + 1;  // number of neighbors to calculate thetas for
	int n_sum, nsteps;       // number of fourier coefficients to use in calculation of theta
	data >> n_sum;
	
	/* dimension work space */
	ifstreamT& in = Input();
	double tstep, totaltime, looptime;  // timestep
	in >> fN_times;
	in >> tstep;
	totaltime = fN_times * tstep;   // total time of simulation
	
	/* determine correct loop time for theta and time history variables */
	if (totaltime <= 3.0)  // assume that t_crit = 3.0 in normalized time
		looptime = totaltime;
	else
		looptime = 3.0;   // need to normalize this time
	fNumstep_crit = int(3.0/tstep) + 1;	// currently assuming normalized t_crit = 3 
										// need to actually compute using LJ parameters
	/* determine correct number of timesteps to store for theta and history variables */
	if (fN_times+1 <= fNumstep_crit)	// WATCH +1 FACTOR!!!
		nsteps = fN_times+1;
	else
		nsteps = fNumstep_crit;

	iArrayT crit(fNumstep_crit);
	fShift.Dimension(fNumstep_crit-1);
	crit.SetValueToPosition();
	fShift.CopyPart(0, crit, 1, fNumstep_crit-1);
	fThetaTable.Dimension(2*fNcrit+1);	//(n_neighbor)
	
	/* dimension boundary atom history */
	fHistoryTablet.Dimension(fTopatoms.Length());   
	fHistoryTableb.Dimension(fBottomatoms.Length());

	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		fHistoryTableb[i].Dimension(nsteps, 2);  // used to be fN_times+1,2
		fHistoryTablet[i].Dimension(nsteps, 2);
	}

	/* dimension theta table data structure */
	for (int i = 0; i < fThetaTable.Length(); i++)
	        fThetaTable[i].Dimension(nsteps, 2*2);  // used to be fN_times+1,2*2
	
	/* read in data table */
	dArrayT row;
	ArrayT<dArray2DT> data_table(n_neighbor);
	for (int i = 0; i < data_table.Length(); i++)
	{
		/* dimension */
		dArray2DT& n_table = data_table[i];
		n_table.Dimension(n_sum, 2*2);
			
		/* read */
		for (int j = 0; j < n_table.MajorDim(); j++)
		{
			int junk;
			data >> junk >> junk;
			
			/* each row is b^T */
			n_table.RowAlias(j, row);
			data >> row;
		}
	}

	dMatrixT theta1(2), temptheta(2), transform(2);
	transform.Identity(1.0);
	transform(1,1) = -1.0;

	/* compute theta's - in correct order for bottom layer, need to transpose for top layer */
	for (int i = 0; i < n_neighbor; i++)	// fThetaTable.Length()
	{
		int count = 0;
		/* get each set of Fourier coefficients */
		const dArray2DT& theta_i = data_table[i];
	       	for (double j = 0.0; j <= looptime; j+=tstep)	// theta only goes to normalized time = 3
		{
			int n_theta = data_table[i].MajorDim();
			temptheta = 0.0;
			for (int k = 0; k < n_theta; k++)  // can truncate this summation
			{
				/* Extract each row of coefficients */
				temptheta.AddScaled(sin((k+1)*pi*j/3.2), theta_i(k)); // tmax = 3.2 now
		      	}

			/* add temptheta into fThetaTable */
			fThetaTable[i].SetRow(count, temptheta);
		
			/* add symmetric portion onto fThetaTable */
			if (i < n_neighbor - 1)
			{
				theta1.MultABC(transform, temptheta, transform);
				fThetaTable[2*fNcrit-i].SetRow(count, theta1);
			}	
			count++;
		}
	}
}

#endif /* BRIDGING_ELEMENT */
