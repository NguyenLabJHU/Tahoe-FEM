/* $Id: FEManagerT_THK.cpp,v 1.2 2007-08-10 18:14:14 d-farrell2 Exp $ */

#include "FEManagerT_THK.h"
#if defined(BRIDGING_ELEMENT)

#include "ifstreamT.h"
#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "TimeManagerT.h"
#include "FieldT.h"
#include "StringT.h"
#include "ParticlePairT.h"
#include "RaggedArray2DT.h"
#include "iGridManagerT.h"
#include "iAutoArrayT.h"
#include "iNodeT.h"
#include "ParameterUtils.h"
#include "ParameterContainerT.h"

#include <iostream.h>
#include <fstream.h>
#include <math.h>

using namespace Tahoe;

const double tol = 1.0e-3;   // for neighbor searching tolerance
const double root32 = sqrt(3.0)/2.0;    // for neighbor searching tolerance

/* constructor */
FEManagerT_THK::FEManagerT_THK(const StringT& input, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	FEManagerT_bridging(input, output, comm, argv, task)
{
	SetName("tahoe_THK");
}

// 2D/3D MD/THK and BSM THK Initialization
void FEManagerT_THK::InitializeTHK(bool ignore_continuum)
{
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
	
	// read other parameters and initialize data
	if (nsd == 2)
		fNeighbors = 2 * fNcrit + 1;	// maximum number of neighbors per atom in 2D (square lattice)
	else if (nsd == 3)
		fNeighbors = (2 * fNcrit + 1)*(2 * fNcrit + 1);	// maximum number of neighbors per atom in 3D (simple cubic)
	else
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize", "%d dimensions not implemented", nsd);

	/* read node set indexes */
	fnumsets = fTHKNodes.Length();				// number of MD THK boundary node sets
	int numsets_temp = fTHKGhostNodes.Length();	// number of MD THK ghost node sets
	
	if (fnumsets != numsets_temp)
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize - # of THK BC sets does not match # of THK Ghost atom Sets");
	
	// collect sets, set up some arrays - here fbound_set_atoms is to be the ghost atoms
	if (fnumsets >= 1)
	{
		fbound_set_atoms.Dimension(fnumsets);
		fghost_set_atoms.Dimension(fnumsets);
		fbound_neighbor_atoms.Dimension(fnumsets);
		for (int i = 0; i < fnumsets; i++)
		{
			// collect sets - minimum 1 - fbound_set_atoms is the ghost plane
			// put each set into array of integer arrays
			fbound_set_atoms[i] = model->NodeSet(fTHKNodes[i]);
			fghost_set_atoms[i] = model->NodeSet(fTHKGhostNodes[i]);
			
			// dimension array of 2D arrays to hold boundary plane neighbors (not in the ghost plane)
			// set the value to -1, for no neighbor
			fbound_neighbor_atoms[i].Dimension(fghost_set_atoms[i].Length(), fNeighbors);
			fbound_neighbor_atoms[i] = -1;
		}
	}
	else
	{
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize - # of THK BC sets less than 1"); // this is the only check on fnumsets
	}
	
	// figure out how many ghost and boundary atoms there are (probably could just ask the nodemanager)
	ftotal_b_atoms = 0;
	ftotal_g_atoms = 0;
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		ftotal_g_atoms += fghost_set_atoms[isetcount].Length();
		ftotal_b_atoms += fbound_set_atoms[isetcount].Length(); 
	}
	
	// dimension array for THK force (sized for 1 entry per ghost atom, no repeats)
	fTHKdisp.Dimension(ftotal_g_atoms, nsd);

	// do neighbor search - only done once as of original implmententation -- by assumptions of method, only needs to be done once
	// this looks for neighboring atoms in the boundary plane (_not_ the ghost plane)
	if (nsd == 2)
		DoNeighSearch2D();
	else if (nsd == 3)
		DoNeighSearch3D();
	else
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize", "%d dimensions not implemented", nsd);
	
	// if needed, obtain the ghost atom properties properties map (to turn off the atoms when computing the FEM force)
	if (ignore_continuum == false) DoGhostMap();
	
	
	if (fTHK_type == "beta")
		ComputeBetaTables();	// compute beta tables
	else if (fTHK_type == "theta")
		ComputeThetaTables();	// compute theta tables
	else
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize", "specified THK type does not match 'beta' or 'theta'");
}

/* return iArrayT of boundary and ghost atom numbers */
const iArrayT& FEManagerT_THK::InterpolationNodes(void)
{
	const iArrayT& ghost = GhostNodes();
	fInterpolationNodes.Dimension(ftotal_g_atoms+ftotal_b_atoms);
	fInterpolationNodes.CopyIn(0,ghost);	// copy in ghost atoms first
	
	
	// loop over each boundary, put boundary atoms into the array(make it so it has no repeats)
	
	int offset = 0;
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		fInterpolationNodes.CopyIn((ghost.Length() + offset), fbound_set_atoms[isetcount]);	// then copy in set atoms
		offset += fbound_set_atoms[isetcount].Length();
	}		

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

//  calculate THK displacement for ghost atoms for 2/3D disp formulation
const dArray2DT& FEManagerT_THK::THKDisp(const StringT& bridging_field, const dArray2DT& badisp)
{
	/* This is the displacement formulation of the THK BC. Here is a rundown of it
	 * The equation is basically u_g = ubar_g + Theta * (u_b - ubar_b)
	 *
	 * 1)	badisp is the coarse scale displacement ubar for the boundary atoms in a
	 *		format with each set sequentially in vector form.
	 * 2)	There are no 'special' atoms. The ghost atom displacement is specified 
	 *		based on the THK and the neighboring boundary atoms so there are no special
	 *		considerations for the corner atoms.
	 * 3)	For now, it is set up for a 'fixed' boundary(for MD/THK - BSM goes as normal).
	 *		Once running, put in second THK and Linear gradient model choices, maybe others.
	 *		all need to do is get badisp somehow, and it will work. Now, badisp is zero
	 */ 
	
	// Get the number of spatial dimensions
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
	
	// dimension some local arrays (these get reused and re defined constantly here)
	dArrayT atomdisp, femdisp(nsd), diff(nsd);
	
	// get some information needed to access the actual MD displacements
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);

	const int stepnum = FEManagerT::StepNumber();	// to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();	// delta t_md

	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, nsd);	// copy all rows of history except last row 

	// Calculate u - ubar for each set
	
	int dispcounter = -1; // keeps track of place in badisp array
	int pos_temp = -1;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// now go through the boundary atom set
		for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
		{
			dispcounter++;
			pos_temp = dispcounter;	// increment the counter
						
			if (stepnum < fNumstep_crit)	// non-shift case
			{
				// for the boundary plane, get the difference q-ubar
				// access the MD displacements directly to save memory
				((*atomfield)[0]).RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
				badisp.RowAlias(pos_temp, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				fHistoryTable[isetcount][i].SetRow(stepnum, diff);
			}
			else	// t > t_crit
			{
				// for the boundary plane
				// access the MD displacements directly to save memory
				((*atomfield)[0]).RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
				badisp.RowAlias(pos_temp, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				shift.RowCollect(fShift, fHistoryTable[isetcount][i]);
				fHistoryTable[isetcount][i].BlockRowCopyAt(shift, 0);
				fHistoryTable[isetcount][i].SetRow(fNumstep_crit-1, diff);
			}
		}
	}
	
	dMatrixT theta(nsd);
	dArrayT gdisp_temp(nsd), disp_temp1(nsd), disp_temp2(nsd);
	InverseMapT setnodes;
	int counter2 = -1;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// set inverse maps  - what does this do?? 
		setnodes.SetMap(fbound_set_atoms[isetcount]);	// Set global to local map
		int dex;

		// calculate the ghost atom fine scale displacement u'_g
		for (int i = 0; i < fbound_neighbor_atoms[isetcount].MajorDim(); i++)	// loop over the ghost node neighbor sets
		{
			counter2++;			// increment the counter
			gdisp_temp = 0.0;	// zero the displacement
			
			for (int count = 0; count < fNeighbors; count++)		// Go through all the neighbors
			{
				if (fbound_neighbor_atoms[isetcount](i,count) != -1)	// if there is a neighbor
				{
					// Get the corresponding displacement history and THK matrix
					dex = setnodes.Map(fbound_neighbor_atoms[isetcount](i,count));
					const dArray2DT& disp_temp0 = fHistoryTable[isetcount][dex];
					const dArray2DT& theta_temp = fThetaTable_array[isetcount][fNeighbors-1-count];
					
					// calculate fine scale THK disp using theta here
					if (stepnum < fNumstep_crit)
					{
						for (int l = 0; l < stepnum; l++)
						{												
							theta.Alias(nsd,nsd,theta_temp(l));
							disp_temp0.RowAlias(stepnum-l, disp_temp1);
							theta.Multx(disp_temp1, disp_temp2);
							gdisp_temp.AddScaled(timestep, disp_temp2);	
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{												
							theta.Alias(nsd,nsd,theta_temp(l));
							disp_temp0.RowAlias(fNumstep_crit-l-1, disp_temp1);
							theta.Multx(disp_temp1, disp_temp2);
							gdisp_temp.AddScaled(timestep, disp_temp2);	
						}
					}
				}
			}
			
			// Set the displacement
			fTHKdisp.SetRow(counter2, gdisp_temp);
		}
	}
	//cout << "fTHKdisp = " << fTHKdisp << endl;
	return fTHKdisp;
}

//  calculate THK displacement for ghost atoms for 2/3D disp formulation, use Beta form
const dArray2DT& FEManagerT_THK::BetaTHKDisp(const StringT& bridging_field, const dArray2DT& badisp, const dArray2DT& bavel)
{
	/* This is the displacement formulation of the Beta THK BC. Here is a rundown of it
	 * The equation is basically u_g = ubar_g + Beta(0)(u_b - ubar_b) - Beta * (v_b - vbar_b)
	 *
	 * 1)	badisp is the coarse scale displacement ubar for the boundary atoms in a
	 *		format with each set sequentially in vector form. bavel is the velocity,
	 *		in the same form.
	 * 2)	There are no 'special' atoms. The ghost atom displacement is specified 
	 *		based on the THK and the neighboring boundary atoms so there are no special
	 *		considerations for the corner atoms.
	 * 3)	For now, it is set up for a 'fixed' boundary(for MD/THK - BSM goes as normal).
	 *		Once running, put in second THK and Linear gradient model choices, maybe others.
	 *		all need to do is get badisp somehow, and it will work. Now, badisp is zero
	 */ 
	
	// Get the number of spatial dimensions
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
	
	// dimension some local arrays (these get reused and re defined constantly here)
	dArrayT atomdisp, atomvel, femdisp(nsd), femvel(nsd), diff(nsd);
	
	// get some information needed to access the actual MD displacements
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);

	const int stepnum = FEManagerT::StepNumber();	// to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();	// delta t_md

	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, nsd);	// copy all rows of history except last row 

	// Calculate (v - vbar) for each set at the current time
	
	int dispcounter = -1; // keeps track of place in bavel array
	int pos_temp = -1;
	
	// declare local array to hold uprime information for all boundary sets at current step.
	// dimension all the way down during loop later.
	ArrayT< ArrayT<dArrayT> > uprime_b;
	uprime_b.Dimension(fnumsets);
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		uprime_b[isetcount].Dimension(fbound_set_atoms[isetcount].Length());
		// now go through the boundary atom set
		for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
		{
			uprime_b[isetcount][i].Dimension(nsd);
			dispcounter++;
			pos_temp = dispcounter;	// increment the counter
						
			if (stepnum < fNumstep_crit)	// non-shift case
			{
				// for the boundary plane, get the difference q-ubar
				// access the MD velocities directly to save memory
				((*atomfield)[1]).RowAlias(fbound_set_atoms[isetcount][i], atomvel);
				bavel.RowAlias(pos_temp, femvel);
				diff.DiffOf(atomvel, femvel);
				fHistoryTable[isetcount][i].SetRow(stepnum, diff);
			}
			else	// t > t_crit
			{
				// for the boundary plane
				// access the MD velocities directly to save memory
				((*atomfield)[1]).RowAlias(fbound_set_atoms[isetcount][i], atomvel);
				bavel.RowAlias(pos_temp, femvel);
				diff.DiffOf(atomvel, femvel);
				shift.RowCollect(fShift, fHistoryTable[isetcount][i]);
				fHistoryTable[isetcount][i].BlockRowCopyAt(shift, 0);
				fHistoryTable[isetcount][i].SetRow(fNumstep_crit-1, diff);
			}
			
			// do the same as above but for displacement, only keep for current time.
			((*atomfield)[0]).RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
			badisp.RowAlias(pos_temp, femdisp);
			uprime_b[isetcount][i].DiffOf(atomdisp, femdisp);
		}
	}
	
	dMatrixT beta(nsd);
	dArrayT gdisp_temp(nsd), disp_temp1(nsd), disp_temp2(nsd), vel_temp1(nsd);
	InverseMapT setnodes;
	int counter2 = -1;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// set inverse maps  - what does this do?? 
		setnodes.SetMap(fbound_set_atoms[isetcount]);	// Set global to local map
		int dex;

		// calculate the ghost atom fine scale displacement u'_g
		for (int i = 0; i < fbound_neighbor_atoms[isetcount].MajorDim(); i++)	// loop over the ghost node neighbor sets
		{
			counter2++;			// increment the counter
			gdisp_temp = 0.0;	// zero the displacement
			
			for (int count = 0; count < fNeighbors; count++)		// Go through all the neighbors
			{
				if (fbound_neighbor_atoms[isetcount](i,count) != -1)	// if there is a neighbor
				{
					// Get the corresponding velocity history, THK matrix and displacement
					dex = setnodes.Map(fbound_neighbor_atoms[isetcount](i,count));
					const dArray2DT& vel_temp0 = fHistoryTable[isetcount][dex];
					const dArray2DT& beta_temp = fThetaTable_array[isetcount][fNeighbors-1-count];
					
					// calculate the first part of the fine scale here (Beta(0)(u_b - ubar_b))
					beta.Alias(nsd,nsd,beta_temp(0));
					beta.Multx(uprime_b[isetcount][dex], disp_temp2);
					gdisp_temp+= disp_temp2;
					
					// calculate fine scale THK disp using beta here (- Beta * (v_b - vbar_b)) and add it on to the above
					if (stepnum < fNumstep_crit)
					{
						for (int l = 0; l <= stepnum; l++)
						{												
							beta.Alias(nsd,nsd,beta_temp(l));
							vel_temp0.RowAlias(stepnum-l, vel_temp1);
							beta.Multx(vel_temp1, disp_temp2);
							gdisp_temp.AddScaled(-timestep, disp_temp2);	
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{												
							beta.Alias(nsd,nsd,beta_temp(l));
							vel_temp0.RowAlias(fNumstep_crit-l-1, vel_temp1);
							beta.Multx(vel_temp1, disp_temp2);
							gdisp_temp.AddScaled(-timestep, disp_temp2);	
						}
					}
				}
			}
			
			// Set the displacement
			fTHKdisp.SetRow(counter2, gdisp_temp);
		}
	}
	//cout << "fTHKdisp = " << fTHKdisp << endl;	
	return fTHKdisp;
}

/* describe the parameters needed by the interface */
void FEManagerT_THK::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FEManagerT_bridging::DefineParameters(list);

	/* time-history kernel parameters */
	list.AddParameter(ParameterT::Integer, "N_crit");
	list.AddParameter(ParameterT::Double, "T_cut");	// user can set cutoff time as desired (maximum is Tmax for the series)
	list.AddParameter(ParameterT::Double, "lattice_parameter");
	list.AddParameter(ParameterT::Double, "interplanar_parameter");	// for finding boundary plane neighbors
	
	ParameterT omega_sys(ParameterT::Double, "Omega_sys");
	omega_sys.SetDefault(1.0);
	list.AddParameter(omega_sys);	// parameter to scale Tcrit and the bn's if the k used in theta file differs
	
	// the ghost mapping for the coupling matrix specified in a file
	list.AddParameter(ParameterT::Word, "ghostmap_file");
	
	// parameter specifying to use theta or beta. default, theta 
	ParameterT THK_type(ParameterT::Word, "THK_type");
	THK_type.SetDefault("theta");
	list.AddParameter(THK_type);
}

/* information about subordinate parameter lists */
void FEManagerT_THK::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FEManagerT_bridging::DefineSubs(sub_list);
	
	// files containing fourier coeffs for THK
	sub_list.AddSub("theta_file_ID_list");

	// nodes affected by THK boundary conditions (boundary atoms)
	sub_list.AddSub("THK_nodes_ID_list");
	
	// nodes affected by THK boundary conditions (ghost atoms)
	sub_list.AddSub("THK_ghost_nodes_ID_list");
	
	// list of outward normals for THK BC planes
	sub_list.AddSub("THK_plane_normals_list", ParameterListT::OnePlus);
	
}

ParameterInterfaceT* FEManagerT_THK::NewSub(const StringT& name) const
{
	if (name == "THK_plane_normals_list")
	{
		ParameterContainerT* x_choice = new ParameterContainerT(name);
		
		/* by dimension */
		x_choice->SetListOrder(ParameterListT::Choice);
		x_choice->AddSub("Vector_2");
		x_choice->AddSub("Vector_3");
	
		return x_choice;
	}
	else // inherited
		return FEManagerT_bridging::NewSub(name);
}

/* accept parameter list */
void FEManagerT_THK::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FEManagerT_bridging::TakeParameterList(list);
	
	int nsd = NodeManager()->NumSD(); // get # spatial dimensions
	
	/* extract THK parameters */
	fNcrit = list.GetParameter("N_crit");
	fTcut = list.GetParameter("T_cut");
	fOmega_sys = list.GetParameter("Omega_sys");
	fLatticeParameter = list.GetParameter("lattice_parameter");
	fSearchParameter = list.GetParameter("interplanar_parameter");
	StringT path;
	path.FilePath(InputFile());
	
	// extract matrix for use when doing ghostoffmap in BridgingScaleManagerT
	fGhostMapFile = list.GetParameter("ghostmap_file");
	fGhostMapFile.ToNativePathName();
	StringT filepath;
	filepath.FilePath(InputFile());
	fGhostMapFile.Prepend(path);
	
	// get THK type - theta or beta
	fTHK_type = list.GetParameter("THK_type");
	
	// Files of fourier coefficients
	const ParameterListT& file_list = list.GetList("theta_file_ID_list");
	StringListT::Extract(file_list, fThetaFile_array);
	for (int i = 0; i < fThetaFile_array.Length(); i++)
	{
		fThetaFile_array[i].ToNativePathName();
		fThetaFile_array[i].Prepend(path);
	}
	

	// nodes affected by THK boundary conditions (boundary atoms)
	const ParameterListT& id_list = list.GetList("THK_nodes_ID_list");
	StringListT::Extract(id_list, fTHKNodes);
	
	// nodes affected by THK boundary conditions (ghost atoms)
	const ParameterListT& ghost_id_list = list.GetList("THK_ghost_nodes_ID_list");
	StringListT::Extract(ghost_id_list, fTHKGhostNodes);
	
	// list of outward normals for THK BC planes
	const ArrayT<ParameterListT>& subs = list.Lists();
	int listcount = 0;
	fTHK_normals.Dimension(fTHKNodes.Length());
	
	for (int i = 0; i < subs.Length(); i++)
	{
		const StringT& name = subs[i].Name();
		if (name == "THK_plane_normals_list" )
		{
			listcount++;
			const ParameterListT& THK_normals_list = subs[i].GetListChoice(*this,"THK_plane_normals_list");
			VectorParameterT::Extract(THK_normals_list, fTHK_normals[listcount-1]);
			if (fTHK_normals[listcount-1].Length() != nsd) 
				ExceptionT::GeneralFail("FEManagerT_THK::TakeParameterList", "\"THK_plane_normals_list\" entry %d should be length %d not %d", i, nsd, fTHK_normals[listcount-1].Length());
		}
		
	}
	if (listcount != fTHKNodes.Length()) 
		ExceptionT::GeneralFail("FEManagerT_THK::TakeParameterList", "\"THK_plane_normals_list\" should be length %d not %d", fTHKNodes.Length(), listcount);
	
}

/*************************************************************************
 * Private
 *************************************************************************/
// perform neighbor search for THK boundary atoms, 2D 
void FEManagerT_THK::DoNeighSearch2D(void)
{
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
	
	// do neighbor search - only done once as of original implmententation -- by assumptions of method, only needs to be done once
	// this looks for neighboring atoms in the boundary plane (_not_ the ghost plane)
	NodeManagerT* node = FEManagerT::NodeManager();				// get the node manager
	const dArray2DT& initcoords = node->InitialCoordinates();	// get initial coords for all atoms
	iArrayT temp_atom(initcoords.MajorDim());					// set up another temp integer array with length = # atoms
	temp_atom.SetValueToPosition();								// initialize the above array to a local numbering scheme
	
	/* configure search grid - CURRENTLY SEARCHING ONLY NON-IMAGE ATOMS (REAL+GHOST) */
	iGridManagerT grid(10, 100, initcoords, &temp_atom);		// set up search grid based on initial coordinates of real + ghost atoms
	grid.Reset();
	dArrayT acoord1, acoord2, ncoord1, ncoord2;					// declare some more temp arrays (maybe we can kill some of these off, save memory?)
	double checkdot, mag1, xmag, ymag;							// declare some temp variables
	int counter;
	dArrayT facenormal;											// nodeset normal (error = 0)											
	facenormal = 0.0;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{ 
		
		// get the normal of plane of boundary atoms (provided in input file)
		facenormal = fTHK_normals[isetcount];
		if (facenormal == 0.0)
			ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D", "THK Boundary face normal not defined");
		
		// find neighbor atoms, based on results of search (2 level search, coarse and fine)
		for (int i = 0; i < fghost_set_atoms[isetcount].Length(); i++)
		{
			initcoords.RowAlias(fghost_set_atoms[isetcount][i], acoord1);		// alias the coordinates of the ith real boundary atom to acoord1
			
			// access the candidate neighbor points
			const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord1.Pointer(), (fNcrit+1)*1.05*sqrt(2.0)*fLatticeParameter); // will work for FCC, BCC, SC 
			
			for (int j = 0; j < hitsa.Length(); j++)		// loop over the candidate neighbors
			{			
				// distance between atom and candidate neighbor
				ncoord1.Alias(nsd, hitsa[j].Coords());		// alias the coordinates of the candidate neighbor to ncoord1
				xmag = ncoord1[0]-acoord1[0];				// get the distance between the atom and candidate neighbor point in each direction
				ymag = ncoord1[1]-acoord1[1];
				mag1 = sqrt(xmag*xmag+ymag*ymag);			// calculate the point to point distance between the neigbor and atom
				
				// now do neighbor search first determine which of the candidates is not in the same plane
				// once the atoms are determined not to be in the same plane, look at the distance, see if it is within the desired range.
				// this will work for arbitrary plane orientations - but lattice parameter is the 'equivalent' square lattice parameter
				// so in the xml file, if a square lattice is being considered, the true lattice parameter is used
				// if a hexagonal lattice is used, one must supply 1/2 the true lattice parameter. 
				// The search is based on square lattice with 'non-atom' points where there is no atom in the case of the hexagonal
							
				// Get in-plane neighbors
				if (mag1 > 1.0e-8) // if mag1 is close enough to be zero -> this would blow up -> but should be zero (allows self checking)
					checkdot = (1/mag1) *((facenormal[0] * xmag) + (facenormal[1] * ymag)); // compute dot product of normal and candidate vector
				else
					checkdot = 0.0;
				
				if (fabs(checkdot) > tol && fabs(mag1*checkdot) < fSearchParameter)	// if the candidate neighbor is not in the ghost plane and within search area, check it
				{																					 					
					// do fine check, store matches in array for THK BC application
					
					// fNCrit denotes neighborshell in the plane
					counter = 0;
					// calculate the second lattice vector (perpendicular to normal) (based on outward normal)
					dArrayT lattice_vect(2);
					lattice_vect[0] = -facenormal[1];
					lattice_vect[1] = facenormal[0];
					
					// find projection of vector to candidate neighbor onto second lattice vector
					double checkdot2 = ((lattice_vect[0] * xmag) + (lattice_vect[1] * ymag));
					
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						// check candidate neighbor position against position in neighbor shell layer
						double dl = l;
						
						if (fabs(checkdot2-(dl*fLatticeParameter)) < tol)
						{																				
							fbound_neighbor_atoms[isetcount](i,counter) = hitsa[j].Tag();									
						}
						counter++;	// counter increments to make sure atoms are in right order																											
					}																				
				}																					
			}																					
		}
	}
}

// perform neighbor search for THK boundary atoms, 3D 
void FEManagerT_THK::DoNeighSearch3D(void)
{
#pragma message("This search is still not very general... figure out better way")
	// This neighbor search curently assumes that the user specifies a normal along one of the coordinate directions
	// So, based on the 2D analogy, it will set up the appropriate triad of lattice vectors to ensure the atoms are 
	// found correctly (basically, assumes that _only_ one of the entries in the normal vector is non-zero)
	
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
	
	// do neighbor search - only done once as of original implmententation -- by assumptions of method, only needs to be done once
	// this looks for neighboring atoms in the boundary plane (_not_ the ghost plane)
	NodeManagerT* node = FEManagerT::NodeManager();				// get the node manager
	const dArray2DT& initcoords = node->InitialCoordinates();	// get initial coords for all atoms
	iArrayT temp_atom(initcoords.MajorDim());					// set up another temp integer array with length = # atoms
	temp_atom.SetValueToPosition();								// initialize the above array to a local numbering scheme
	
	/* configure search grid - CURRENTLY SEARCHING ONLY NON-IMAGE ATOMS (REAL+GHOST) */
	iGridManagerT grid(10, 100, initcoords, &temp_atom);		// set up search grid based on initial coordinates of real + ghost atoms
	grid.Reset();
	dArrayT acoord1, acoord2, ncoord1, ncoord2;					// declare some more temp arrays (maybe we can kill some of these off, save memory?)
	double checkdot, mag1, xmag, ymag, zmag;					// declare some temp variables
	int counter;
	dArrayT facenormal, lat_vect1(3), lat_vect2(3);					// nodeset normal (error = 0) and other lattice vectors											
	facenormal = 0.0;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{ 
		
		// get the normal of plane of boundary atoms (provided in input file)
		facenormal = fTHK_normals[isetcount];
		if (facenormal == 0.0)
			ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D", "THK Boundary face normal not defined");
		
		if (facenormal[0] == 1.0)	// outward normal in positive x-direction
		{
			lat_vect1[0] = 0.0;		// lattice vector in y-direction
			lat_vect1[1] = -1.0;
			lat_vect1[2] = 0.0;
			
			lat_vect2[0] = 0.0;		// lattice vector in z-direction
			lat_vect2[1] = 0.0;
			lat_vect2[2] = 1.0; 
		}
		else if (facenormal[0] == -1.0)	// outward normal in negative x-direction
		{
			lat_vect1[0] = 0.0;		// lattice vector in y-direction
			lat_vect1[1] = 1.0;
			lat_vect1[2] = 0.0;
			
			lat_vect2[0] = 0.0;		// lattice vector in z-direction
			lat_vect2[1] = 0.0;
			lat_vect2[2] = 1.0; 
		}
		else if (facenormal[1] == 1.0)	// outward normal in positive y-direction
		{
			lat_vect1[0] = 1.0;		// lattice vector in x-direction
			lat_vect1[1] = 0.0;
			lat_vect1[2] = 0.0;
			
			lat_vect2[0] = 0.0;		// lattice vector in z-direction
			lat_vect2[1] = 0.0;
			lat_vect2[2] = 1.0; 
		}
		else if (facenormal[1] == -1.0)	// outward normal in negative y-direction
		{
			lat_vect1[0] = -1.0;		// lattice vector in x-direction
			lat_vect1[1] = 0.0;
			lat_vect1[2] = 0.0;
			
			lat_vect2[0] = 0.0;		// lattice vector in z-direction
			lat_vect2[1] = 0.0;
			lat_vect2[2] = 1.0; 
		}
		else if (facenormal[2] == 1.0)	// outward normal in positive z-direction
		{
			lat_vect1[0] = 0.0;		// lattice vector in y-direction
			lat_vect1[1] = -1.0;
			lat_vect1[2] = 0.0;
			
			lat_vect2[0] = 1.0;		// lattice vector in x-direction
			lat_vect2[1] = 0.0;
			lat_vect2[2] = 0.0; 
		}
		else if (facenormal[2] == -1.0)	// outward normal in negative z-direction
		{
			lat_vect1[0] = 0.0;		// lattice vector in y-direction
			lat_vect1[1] = 1.0;
			lat_vect1[2] = 0.0;
			
			lat_vect2[0] = 1.0;		// lattice vector in x-direction
			lat_vect2[1] = 0.0;
			lat_vect2[2] = 0.0; 
		}
		 		
		// find in-plane neighbor atoms, based on results of search (2 level search, coarse and fine)
		for (int i = 0; i < fghost_set_atoms[isetcount].Length(); i++)
		{
			initcoords.RowAlias(fghost_set_atoms[isetcount][i], acoord1);		// alias the coordinates of the ith real boundary atom to acoord1
			
			// access the candidate neighbor points
			const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord1.Pointer(), (fNcrit+1)*1.05*sqrt(2.0)*fLatticeParameter); // will work for FCC, BCC, SC 
																																// DEF added (fNCrit+1) to scale the search area
			for (int j = 0; j < hitsa.Length(); j++)		// loop over the candidate neighbors
			{			
				// distance between atom and candidate neighbor
				ncoord1.Alias(nsd, hitsa[j].Coords());		// alias the coordinates of the candidate neighbor to ncoord1 (perhaps change syntax?)
				xmag = ncoord1[0]-acoord1[0];				// get the distance between the atom and candidate neighbor point in each direction
				ymag = ncoord1[1]-acoord1[1];
				zmag = ncoord1[2]-acoord1[2];
				mag1 = sqrt(xmag*xmag+ymag*ymag+zmag*zmag);	// calculate the point to point distance between the neigbor and atom
				
				// now do neighbor search first determine which of the candidates is in the same plane (look to see if facenormal & normalized vector between atoms is ~ +/- tol)
				// once the atoms are determined to be in the same plane, look at the distance, see if it is within the desired range.
				// this will work for arbitrary plane orientations - but lattice parameter is the true parameter
							
				// Get boundary plane neighbors
				if (mag1 > 1.0e-8) // if mag1 is close enough to be zero -> this would blow up -> but should be zero (allows self checking)
					checkdot = (1/mag1) *((facenormal[0] * xmag) + (facenormal[1] * ymag) + (facenormal[2] * zmag)); // compute dot product of normal and candidate vector
				else
					checkdot = 0.0;
				
				if (fabs(checkdot) > tol && fabs(mag1*checkdot) < fSearchParameter)	// if the candidate neighbor is not in the ghost plane, and within search area check it
				{																					 					
					// do fine check, store matches in array for THK BC application
					// Project the vector to the candidate onto the lattice vectors from above
					double checkdot1 = ((lat_vect1[0] * xmag) + (lat_vect1[1] * ymag) + (lat_vect1[2] * zmag));
					double checkdot2 = ((lat_vect2[0] * xmag) + (lat_vect2[1] * ymag) + (lat_vect2[2] * zmag));
					
					// fNCrit denotes neighborshell in the boundary plane
					counter = 0;
					for (int k = -fNcrit; k <= fNcrit; k++)
					{
						double dk = k;
						for (int l = -fNcrit; l<= fNcrit; l++)
						{
							double dl = l;
							if ( fabs(checkdot1-(dk*fLatticeParameter)) < tol && fabs(checkdot2-(dl*fLatticeParameter)) < tol)	// check position in neighbor shell, assumes cubic
							{																				
								fbound_neighbor_atoms[isetcount](i,counter) = hitsa[j].Tag();									
							}
							counter++;		// counter increments to make sure atoms are in right order
						}																			
					}																				
				}																					
			}																					
		}
	}
}

// find the ghost atom properties map
void FEManagerT_THK::DoGhostMap(void)
{
	// read in the ghostoffmap into fghostoffmap
	ifstreamT data2(fGhostMapFile);
	if (!data2.is_open())
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D", "file not found: %s", fGhostMapFile.Pointer());
	
	// get dimension of the file from the properties map
	nMatrixT<int>& promap = PropertiesMap(0);   // element group for particles = 0
	const int promap_dim = promap.Rows(); // assumes square property matrix
	// ghostoffmap matrix
	fghostoffmap.Dimension(promap_dim);
	fghostoffmap = promap;
	
	int ghosti = 0; // row number
	int ghostj = 0; // column number
	int ghostval=0;	// value for zero-force potential
	// now read in the rows and columns of the fghostoffmap 2D array, and the value to set it to
	/* Here is how this works:
	 * The ghostonmap is the same as the properties map. The properties map is what is specified in the input file. Each line
	 * from the soecified pair_particle_interaction's is a new property, where the first = 0, etc. This is the upper triangle
	 * of the properties map, which is symmetric. For the ghostoffmap, you want to turn off all of the interactions that include
	 * the ghost atoms. To do this, you must set the entry corresponding to the atom type to a interaction value which will produce
	 * no force (i.e. k=0 for harmonic, epsilon = 0 for LJ). the ghostoffmapfile contains the coordinates for the entries which need
	 * to change, as well as the value they need to be set to. all the other interactions remain the same as the input file.
	 */ 
	for (int ghostk = 0; ghostk < 2*(promap_dim-1); ghostk++)
	{
		data2 >> ghosti >> ghostj >> ghostval ;
		fghostoffmap(ghosti,ghostj) = ghostval;
	}
}

// compute theta tables for 2D/3D disp/disp or disp/force formulation (doesn't matter, its all the same)
void FEManagerT_THK::ComputeThetaTables(void)
{	
	const char caller[] = "FEManagerT_THK::ComputeThetaTables";
	
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
		
	/* dimensions */
	double pi = acos(-1.0);
	int nsteps;       // number of steps to be stored
	
	/* dimension work space */
	const TimeManagerT* time_manager = TimeManager();	// get the time manager
	fN_times = time_manager->NumberOfSteps();		 	// get the number of timesteps
	double tstep = time_manager->TimeStep(); 		 	// get the timestep size
	double totaltime = fN_times * tstep; 				// total time of simulation
	
	/* determine correct loop time for theta and time history variables */
	double looptime = 0.0;
	if (totaltime <= fTcut/fOmega_sys)  // need to store/calculate up to this time (time cutoff)
		looptime = totaltime;
	else
		looptime = fTcut/fOmega_sys;   // need to store/calculate up to this time

	fNumstep_crit = int(((fTcut/fOmega_sys)/tstep) + 0.5) + 1;	// determine number of steps to t_cut
	// DEF note: the 0.5 was added on to ensure that the int chop gets the correct number (to avoid roundoff issues)
	
	/* determine correct number of timesteps to store for theta and history variables */
	if (fN_times+1 <= fNumstep_crit)
		nsteps = fN_times+1;  
	else
		nsteps = fNumstep_crit;

	iArrayT crit(fNumstep_crit);
	fShift.Dimension(fNumstep_crit-1);
	crit.SetValueToPosition();
	fShift.CopyPart(0, crit, 1, fNumstep_crit-1);
	
	// Dimension the arrays which will hold some other arrays (dimension to # of BC sets)
	fThetaTable_array.Dimension(fnumsets);
	fHistoryTable.Dimension(fnumsets);
	dArrayT row;
	ArrayT< ArrayT<dArray2DT> > data_table_array(fnumsets);
	
	int junki;	// trash variable -> dummy for reading in THK data file
	StringT junks;
	
	int n_sum;	// number of modes to use in the THK series (from file)
	double t_max,laplace_a;	// Tmax and a from laplace inverse for THK (from file)
		
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// Open the file with the fourier coeffs
		const StringT& data_file = fThetaFile_array[isetcount];
		ifstreamT data(data_file);
		if (!data.is_open())
			ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
		// crude file reading to get data
		data >> junks >> junks >> junks >> junks >> junks >> junks >> junks >> t_max >> junks >> n_sum >> junks >> laplace_a;	// read in the data needed
		data >> junks >> junks >> junks >> junks >> junks >> junks;	// skip the header
		
		if (t_max < fTcut)		// if the specified cutoff is too long - error out
			ExceptionT::GeneralFail(caller, "T_cut is larger than Tmax: not a valid time cutoff");
		
		// Dimension some arrays of arrays to hold information (dimension to # neighbors, # boundary atoms)
		fThetaTable_array[isetcount].Dimension(fNeighbors);
		fHistoryTable[isetcount].Dimension(fbound_set_atoms[isetcount].Length());
		
		// Now dimension all the way down
		for (int count1 = 0; count1 < fNeighbors; count1++)
		{
			fThetaTable_array[isetcount][count1].Dimension(nsteps, nsd*nsd); // Dimension to hold elements of THK matrix for each of nsteps
		}
		for (int count1 = 0; count1 < fbound_set_atoms[isetcount].Length(); count1++)
		{
			fHistoryTable[isetcount][count1].Dimension(nsteps, nsd); // Dimension to hold displacement history for each of nsteps
		}
		
		// read in data tables of Fourier coefficients for calculation of Theta
		data_table_array[isetcount].Dimension(fNeighbors);	// dimension for each set
		for (int i = 0; i < data_table_array[isetcount].Length(); i++)
		{
			/* dimension */
			dArray2DT& n_table = data_table_array[isetcount][i];
			n_table.Dimension(n_sum, nsd*nsd);
					
			/* read */
			for (int j = 0; j < n_table.MajorDim(); j++)
			{
				if (nsd == 2)
					data >> junki >> junki ; // first 2 cols in file are not needed
				else if (nsd == 3 )
					data >> junki >> junki >> junki ; // first 3 cols in file are not needed
				
				/* each row is b^T */
				n_table.RowAlias(j, row);	// read in the values of the fourier coefficients (coefficients from fourier series expansion)
				data >> row;				// See below for the use of the expansion to calculate theta from the coefficients
			}								// These are simply the fourier sine coeffs for each atom in the unit cell
		}
		
		
		// Now use fourier sine series to compute Theta -> this will end up recycling memory
		dMatrixT theta1(nsd), temptheta(nsd); 

		// compute theta's - uses fourier sine series to represent the quantity
		for (int i = 0; i < fNeighbors; i++)
		{
			int count = 0;
			
			/* get each set of Fourier coefficients */
			const dArray2DT& theta_i = data_table_array[isetcount][i];
			
			for (double j = 0.0; j <= looptime; j+=tstep)	// theta only goes to normalized time = looptime
			{
				int n_theta = data_table_array[isetcount][i].MajorDim();
				temptheta = 0.0;
				
				for (int k = 0; k < n_theta; k++)  // loop over each term in the series
				{
					/* Extract each row of coefficients */
					temptheta.AddScaled(fOmega_sys*exp(laplace_a*j*fOmega_sys)*sin((k+1)*pi*j/(t_max/fOmega_sys)), theta_i(k));
				}

				/* add temptheta into fThetaTable */
				fThetaTable_array[isetcount][i].SetRow(count, temptheta);				
				count++;
			}
			
// This is to use the exact analytical kernel
//			int n_theta = data_table_array[isetcount][i].MajorDim();
//			for (int k = 0; k < n_theta; k++)	// n_theta = n_steps
//			{
//				temptheta = 0.0;
//
//				/* add temptheta into fThetaTable */
//				fThetaTable_array[isetcount][i].SetRow(count, theta_i(count));				
//				count++;
//			}			
		}
	}
}

// compute beta tables for 2D/3D disp/vel or force/vel formulation (doesn't matter, its all the same)
// note - in here, I make use of the same variable names as for the theta case - this is for ease. 
// Make no mistake - this calculates beta. This also uses the same file format as in the case of theta
void FEManagerT_THK::ComputeBetaTables(void)
{	
	const char caller[] = "FEManagerT_THK::ComputeBetaTables";
	
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
		
	/* dimensions */
	double pi = acos(-1.0);
	int nsteps;       // number of steps to be stored
	
	/* dimension work space */
	const TimeManagerT* time_manager = TimeManager();	// get the time manager
	fN_times = time_manager->NumberOfSteps();		 	// get the number of timesteps
	double tstep = time_manager->TimeStep(); 		 	// get the timestep size
	double totaltime = fN_times * tstep; 				// total time of simulation
	
	/* determine correct loop time for theta and time history variables */
	double looptime = 0.0;
	if (totaltime <= fTcut/fOmega_sys)  // need to store/calculate up to this time (time cutoff)
		looptime = totaltime;
	else
		looptime = fTcut/fOmega_sys;   // need to store/calculate up to this time

	fNumstep_crit = int(((fTcut/fOmega_sys)/tstep) + 0.5) + 1;	// determine number of steps to t_cut
	// DEF note: the 0.5 was added on to ensure that the int chop gets the correct number (to avoid roundoff issues)
	
	/* determine correct number of timesteps to store for theta and history variables */
	if (fN_times+1 <= fNumstep_crit)
		nsteps = fN_times+1;  
	else
		nsteps = fNumstep_crit;

	iArrayT crit(fNumstep_crit);
	fShift.Dimension(fNumstep_crit-1);
	crit.SetValueToPosition();
	fShift.CopyPart(0, crit, 1, fNumstep_crit-1);
	
	// Dimension the arrays which will hold some other arrays (dimension to # of BC sets)
	fThetaTable_array.Dimension(fnumsets);
	fHistoryTable.Dimension(fnumsets);
	dArrayT row;
	ArrayT< ArrayT<dArray2DT> > data_table_array(fnumsets);
	
	int junki;	// trash variable -> dummy for reading in THK data file
	StringT junks;
	
	int n_sum;	// number of modes to use in the THK series (from file) - in this case, includes 0 mode.
	double t_max,laplace_a;	// Tmax and a from laplace inverse for THK (from file)
		
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// Open the file with the fourier coeffs
		const StringT& data_file = fThetaFile_array[isetcount];
		ifstreamT data(data_file);
		if (!data.is_open())
			ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
		// crude file reading to get data
		data >> junks >> junks >> junks >> junks >> junks >> junks >> junks >> t_max >> junks >> n_sum >> junks >> laplace_a;	// read in the data needed
		data >> junks >> junks >> junks >> junks >> junks >> junks;	// skip the header
		
		if (t_max < fTcut)		// if the specified cutoff is too long - error out
			ExceptionT::GeneralFail(caller, "T_cut is larger than Tmax: not a valid time cutoff");
		
		// Dimension some arrays of arrays to hold information (dimension to # neighbors, # boundary atoms)
		fThetaTable_array[isetcount].Dimension(fNeighbors);
		fHistoryTable[isetcount].Dimension(fbound_set_atoms[isetcount].Length());
		
		// Now dimension all the way down
		for (int count1 = 0; count1 < fNeighbors; count1++)
		{
			fThetaTable_array[isetcount][count1].Dimension(nsteps, nsd*nsd); // Dimension to hold elements of THK matrix for each of nsteps
		}
		for (int count1 = 0; count1 < fbound_set_atoms[isetcount].Length(); count1++)
		{
			fHistoryTable[isetcount][count1].Dimension(nsteps, nsd); // Dimension to hold velocity history for each of nsteps
		}
		
		// read in data tables of Fourier coefficients for calculation of Theta
		data_table_array[isetcount].Dimension(fNeighbors);	// dimension for each set
		for (int i = 0; i < data_table_array[isetcount].Length(); i++)
		{
			/* dimension */
			dArray2DT& n_table = data_table_array[isetcount][i];
			n_table.Dimension(n_sum, nsd*nsd);
					
			/* read */
			for (int j = 0; j < n_table.MajorDim(); j++)
			{
				if (nsd == 2)
					data >> junki >> junki ; // first 2 cols in file are not needed
				else if (nsd == 3 )
					data >> junki >> junki >> junki ; // first 3 cols in file are not needed
				
				/* each row is b^T */
				n_table.RowAlias(j, row);	// read in the values of the fourier coefficients (coefficients from fourier series expansion)
				data >> row;				// See below for the use of the expansion to calculate beta from the coefficients
			}								// These are simply the fourier cosine coeffs for each atom in the unit cell
		}
		
		
		// Now use fourier sine series to compute Theta -> this will end up recycling memory
		dMatrixT theta1(nsd), temptheta(nsd); 

		// compute beta's - uses fourier cosine series to represent the quantity
		for (int i = 0; i < fNeighbors; i++)
		{
			int count = 0;
			
			/* get each set of Fourier coefficients */
			const dArray2DT& theta_i = data_table_array[isetcount][i];
			
			for (double j = 0.0; j <= looptime; j+=tstep)	// beta only goes to normalized time = looptime
			{
				int n_theta = data_table_array[isetcount][i].MajorDim();
				temptheta = 0.0;
				
				for (int k = 0; k < n_theta; k++) // loop over each term in the series
				{
					/* Extract each row of coefficients */
					temptheta.AddScaled(cos(k*pi*j/(t_max/fOmega_sys)), theta_i(k));
				}

				/* add temptheta into fThetaTable */
				fThetaTable_array[isetcount][i].SetRow(count, temptheta);				
				count++;
			}
			
// This is to use the exact analytical kernel
//			int n_theta = data_table_array[isetcount][i].MajorDim();
//			for (int k = 0; k < n_theta; k++)	// n_theta = n_steps
//			{
//				temptheta = 0.0;
//
//				/* add temptheta into fThetaTable */
//				fThetaTable_array[isetcount][i].SetRow(count, theta_i(count));				
//				count++;
//			}			
		}
	}
}

#endif /* BRIDGING_ELEMENT */
