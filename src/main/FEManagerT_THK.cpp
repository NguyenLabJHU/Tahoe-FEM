/* $Id: FEManagerT_THK.cpp,v 1.21 2005-04-06 15:42:00 paklein Exp $ */
#include "FEManagerT_THK.h"
#if defined(BRIDGING_ELEMENT) && defined(BRIDGING_ELEMENT_DEV)

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

#include <iostream.h>
#include <fstream.h>
#include <math.h>

using namespace Tahoe;

const double tol = 1.0e-3;   // for neighbor searching tolerance 1.0e-8 originally
const double root32 = sqrt(3.0)/2.0;    // for neighbor searching tolerance

/* constructor */
FEManagerT_THK::FEManagerT_THK(const StringT& input, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, TaskT task):
	FEManagerT_bridging(input, output, comm, argv, task)
{
	SetName("tahoe_THK");
}

/* 2D Bridging Scale Initialization */
void FEManagerT_THK::Initialize2D(void)
{
	ModelManagerT* model = FEManagerT::ModelManager();

	// read other parameters and initialize data
	int num_neighbors = 2 * fNcrit + 1;   // maximum number of neighbors per atom in 2D

	/* read node set indexes */
	fnumsets = fTHKNodes.Length();		// number of MD THK boundary node sets
 
	/* collect sets - currently assuming exactly 2 MD THK boundaries in 2D */
	/* further assumes bottom is first node set, top is second node set */
	fBottomatoms = model->NodeSet(fTHKNodes[0]);
	fTopatoms = model->NodeSet(fTHKNodes[1]);
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
	
///// May wish to move this into another method	 NOTE: the file is formatted in rows with 2 columns (i,j) coords in mapping array
	// read in the ghostoffmap into fghostoffmap
	ifstreamT data(fGhostMapFile);
	if (!data.is_open())
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D", "file not found: %s", fGhostMapFile.Pointer());
	
	// get dimension of the file from the properties map
	nMatrixT<int>& promap = PropertiesMap(0);   // element group for particles = 0
	const int promap_dim = promap.Rows(); // assumes square property matrix
	
	int ghosti = 0; // row number
	int ghostj = 0; // column number
	// now read in the rows and columns of the fghostoffmap 2D array
	for (int ghostk = 0; ghostk < promap_dim; ghostk++)
	{
		data >> ghosti >> ghostj ;
		fghostoffmap(ghosti,ghostj) = 1;
	}

///////	
	
	
	/* compute theta tables */
	ComputeThetaTables2D(fThetaFile);
}

/* Bridging Scale 3D Initialization */
void FEManagerT_THK::Initialize3D(void)
{
	/* Implement 3D version of initialize here */
	/* Find neighbors of top 2 planes of atoms */
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();

	// read other parameters and initialize data
	int num_neighborsx = 2 * fNcrit + 1;	// assume same number of neighbors in both x and y directions
	int num_neighborsy = 2 * fNcrit + 1;
	fNeighbors = num_neighborsx * num_neighborsy;

	/* read node set indexes */
	fnumsets = fTHKNodes.Length();		// number of MD THK boundary node sets

/////////////// This whole thing looks like it needs to be generalized
#pragma message("FEManager_THK::Initialize3D only set up for bottom or bottom/top boundary")
 
	// collect sets - currently assuming exactly 1 or 2 MD THK boundaries in 3D
	// further assumes bottom is first node set or only set, top is second node set
	fBottomatoms = model->NodeSet(fTHKNodes[0]);
	fBottomrow.Dimension(fBottomatoms.Length());
	fBottomrow.SetValueToPosition();
	fBottom.Dimension(fBottomatoms.Length(), fNeighbors);
	fBottom = -1;  // -1 -> no neighbor
	
	if (fnumsets == 2)
	{
		fTopatoms = model->NodeSet(fTHKNodes[1]);
		fToprow.Dimension(fTopatoms.Length());
		fToprow.SetValueToPosition();
		fToprow+=fBottomrow.Length();
		fTop.Dimension(fTopatoms.Length(), fNeighbors);
		fTop = -1;  // -1 -> no neighbor
	}
	else if (fnumsets > 2)
	{
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D - only set up for bottom or bottom/top boundaries - tell dave");
	}

	NodeManagerT* node = FEManagerT::NodeManager();
	const dArray2DT& initcoords = node->InitialCoordinates();  // init coords for all atoms
	iArrayT asdf(initcoords.MajorDim());
	asdf.SetValueToPosition();
	
	/* configure search grid - CURRENTLY SEARCHING ONLY NON-IMAGE ATOMS (REAL+GHOST) */
	iGridManagerT grid(10, 100, initcoords, &asdf);
	grid.Reset();
	dArrayT acoord1, acoord2, ncoord1, ncoord2;
	double xdist1, ydist1, mag1, xmag, ymag, zmag;
	int blah2;
	
	// find bottom plane0 neighbor atoms -> at minimum
	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		initcoords.RowAlias(fBottomatoms[i], acoord1);
		/* candidate points for row0 bottom atoms */
		const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord1.Pointer(), 1.05*sqrt(2.0)*fLatticeParameter);
	
		for (int j = 0; j < hitsa.Length(); j++)
		{			
			/* distance between atoms */
			ncoord1.Alias(nsd, hitsa[j].Coords());
			xmag = ncoord1[0]-acoord1[0];
			ymag = ncoord1[1]-acoord1[1];
			zmag = ncoord1[2]-acoord1[2];
			mag1 = sqrt(xmag*xmag+ymag*ymag+zmag*zmag);
						
			/* Get row0 neighbors */
			if (fabs(acoord1[2]-ncoord1[2]) < tol && mag1 < 1.03*sqrt(2.0)*fLatticeParameter)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = fabs(-ncoord1[0]+acoord1[0]+fLatticeParameter*k);
						ydist1 = fabs(-ncoord1[1]+acoord1[1]+fLatticeParameter*l);
						if (xdist1 < tol && ydist1 < tol)
							fBottom(i,blah2) = hitsa[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}
	
	if (fnumsets == 2) // do top plane
	{
		/* find top plane0 neighbor atoms */
		for (int i = 0; i < fTopatoms.Length(); i++)
		{
			initcoords.RowAlias(fTopatoms[i], acoord2);
			/* candidate points for row0 top atoms */
			const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord2.Pointer(), 1.05*sqrt(2.0)*fLatticeParameter);
		
			for (int j = 0; j < hitsa.Length(); j++)
			{			
				/* distance between atoms */
				ncoord2.Alias(nsd, hitsa[j].Coords());
				xmag = ncoord2[0]-acoord2[0];
				ymag = ncoord2[1]-acoord2[1];
				zmag = ncoord2[2]-acoord2[2];
				mag1 = sqrt(xmag*xmag+ymag*ymag+zmag*zmag);
							
				/* Get row0 neighbors */
				if (fabs(acoord2[2]-ncoord2[2]) < tol && mag1 < 1.03*sqrt(2.0)*fLatticeParameter)	
				{
					/* now sort into x by y array for THK BC application */
					blah2 = 0;
					for (int k = -fNcrit; k <= fNcrit; k++)
					{
						for (int l = -fNcrit; l<= fNcrit; l++)
						{
							xdist1 = fabs(-ncoord2[0]+acoord2[0]+fLatticeParameter*k);
							ydist1 = fabs(-ncoord2[1]+acoord2[1]+fLatticeParameter*l);
							if (xdist1 < tol && ydist1 < tol)
								fTop(i,blah2) = hitsa[j].Tag();
							blah2++;
						}	
					}
				}
			}
		}
	}

	
///// May wish to move this into another method	 NOTE: the file is formatted in rows with 2 columns (i,j) coords in mapping array
	// read in the ghostoffmap into fghostoffmap
	ifstreamT data(fGhostMapFile);
	if (!data.is_open())
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D", "file not found: %s", fGhostMapFile.Pointer());
	
	// get dimension of the file from the properties map
	nMatrixT<int>& promap = PropertiesMap(0);   // element group for particles = 0
	const int promap_dim = promap.Rows(); // assumes square property matrix
	// ghostoffmap matrix
	fghostoffmap.Dimension(promap_dim);
	fghostoffmap = 0;
	
	int ghosti = 0; // row number
	int ghostj = 0; // column number
	// now read in the rows and columns of the fghostoffmap 2D array
	for (int ghostk = 0; ghostk < promap_dim; ghostk++)
	{
		data >> ghosti >> ghostj ;
		fghostoffmap(ghosti,ghostj) = 1;
	}

///////	
	
	
	/* Compute Theta tables - may need to pass multiple data files as arguments */
	ComputeThetaTables3D(fThetaFile);
}

/* return iArrayT of boundary and ghost atom numbers - 3D version */
const iArrayT& FEManagerT_THK::InterpolationNodes3D(void)
{
	const iArrayT& ghost = GhostNodes();
	fInterpolationNodes.Dimension(ghost.Length()+fBottomatoms.Length()+fTopatoms.Length());
	fInterpolationNodes.CopyIn(0,ghost);	// copy in ghost atoms first
	fInterpolationNodes.CopyIn(ghost.Length(), fBottomatoms);	// then copy in bottom atoms
	fInterpolationNodes.CopyIn(ghost.Length()+fBottomatoms.Length(), fTopatoms);	// finally top atoms
	return fInterpolationNodes;
}

/* return iArrayT of boundary and ghost atom numbers - 2D version */
const iArrayT& FEManagerT_THK::InterpolationNodes2D(void)
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

/* calculate external force on MD boundary atoms for 2D disp/force formulation */
const dArray2DT& FEManagerT_THK::THKForce(const StringT& bridging_field, const dArray2DT& badisp)
{
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
						theta.Alias(2,2,theta_b(k));
						disp_b.RowAlias(stepnum-k, force1);
						theta.Multx(force1, force2);
						force3b.AddScaled(timestep, force2);
					
						/* top row */
						theta.Alias(2,2,theta_t(k));
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
						theta.Alias(2,2,theta_b(k));
						disp_b.RowAlias(fNumstep_crit-k-1, force1);
						theta.Multx(force1, force2);
						force3b.AddScaled(timestep, force2);
			
						/* top row */
						theta.Alias(2,2,theta_t(k));
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

/* calculate impedance force using 3D disp/force formulation */
const dArray2DT& FEManagerT_THK::THKDisp(const StringT& bridging_field, const dArray2DT& badisp)
{
	/* badisp is in format bottom displacements first, then top */
	/* COULD DIMENSION ONCE IN INTERPOLATIONNODES3D() */
	fTHKforce.Dimension(badisp.MajorDim(), 3);  
	fTHKforce = 0.0;

	/* sort top and bottom displacements into separate arrays - may not be necessary */
	/* ADD topdisp0, etc as class variables so only dimension once? */
	/* COULD DIMENSION ONCE IN INTERPOLATIONNODES3D() */
	dArray2DT topdisp0, bottomdisp0;
	dArrayT atomdisp, femdisp(3), diff(3);
	bottomdisp0.Dimension(fBottomatoms.Length(), 3);
	topdisp0.Dimension(fTopatoms.Length(), 3);
	bottomdisp0.RowCollect(fBottomrow, badisp);
	topdisp0.RowCollect(fToprow, badisp);
	
	/* access the actual MD displacements */
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);
	dArray2DT mddisp = (*atomfield)[0];

	const int stepnum = FEManagerT::StepNumber();  // to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();  // delta t_md

	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, 3);	// copy all rows of history except last row 

	/* Calculate q - ubar for bottom or top/bottom plane 0 of atoms */
	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		/* non-shift case */
		if (stepnum < fNumstep_crit)   
		{
			if (fnumsets == 2)
			{
				// Row 0 top plane of atoms
				mddisp.RowAlias(fTopatoms[i], atomdisp);
				topdisp0.RowAlias(i, femdisp);	
				diff.DiffOf(atomdisp, femdisp);
				fHistoryTablet[i].SetRow(stepnum, diff);
			}
		
			// Row 0 bottom plane of atoms
			mddisp.RowAlias(fBottomatoms[i], atomdisp);
			bottomdisp0.RowAlias(i, femdisp);	
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTableb[i].SetRow(stepnum, diff);
		}
		else	// t > t_crit
		{
			if (fnumsets == 2)
			{
				// Row 0 top plane of atoms
				mddisp.RowAlias(fTopatoms[i], atomdisp);
				topdisp0.RowAlias(i, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				shift.RowCollect(fShift, fHistoryTablet[i]);
				fHistoryTablet[i].BlockRowCopyAt(shift, 0);
				fHistoryTablet[i].SetRow(fNumstep_crit-1, diff);
			}
	
			/* Row 0 bottom plane of atoms */
			mddisp.RowAlias(fBottomatoms[i], atomdisp);
			bottomdisp0.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTableb[i]);
			fHistoryTableb[i].BlockRowCopyAt(shift, 0);
			fHistoryTableb[i].SetRow(fNumstep_crit-1, diff);
		}
	}
		
	dMatrixT theta(3);
	dArrayT force1(3), force2(3), force0a(3), force0b(3);
	
	/* set inverse maps */
	InverseMapT topnodes0, bottomnodes0;
	if (fnumsets == 2) { topnodes0.SetMap(fTopatoms); }	// Set global to local map for top atoms plane 0
	bottomnodes0.SetMap(fBottomatoms);	// Set global to local map for bottom atoms plane 0
	int count, dex;

	/* calculate THK force for top/bottom plane0 atoms */
	for (int i = 0; i < fBottom.MajorDim(); i++)
	{
		force0a = 0.0;
		force0b = 0.0;
		count = 0;
		for (int j = 0; j < 2*fNcrit+1; j++)
		{
			for (int k = 0; k < 2*fNcrit+1; k++)
			{
				if (fBottom(i,count) == -1)	
					int blah = 0;
				else
				{
					/* bottom plane 0 */
					dex = bottomnodes0.Map(fBottom(i,count));
					const dArray2DT& disp_b0 = fHistoryTableb[dex];
					const dArray2DT& theta_b0 = fThetaTableB[fNeighbors-1-count]; 	// start count at 8, go backwards 
					
					/* calculate fine scale THK disp using theta here */
					if (stepnum < fNumstep_crit)  // < fNumstep_crit 
					{
						for (int l = 0; l < stepnum; l++)	// ORIGINALLY l < stepnum
						{												
							/* bottom plane 0 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{												
							/* bottom plane 0 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	
						}
					}
				}
				if (fnumsets == 2)
				{
					if (fTop(i,count) == -1 )	
						int blah = 0;
					else
					{
						/* top plane 0 */
						dex = topnodes0.Map(fTop(i,count));
						const dArray2DT& disp_t0 = fHistoryTablet[dex];
						const dArray2DT& theta_t0 = fThetaTableT[fNeighbors-1-count]; 	// start count at 8, go backwards 
						
						/* calculate fine scale THK disp using theta here */
						if (stepnum < fNumstep_crit)  // < fNumstep_crit 
						{
							for (int l = 0; l < stepnum; l++)	// ORIGINALLY l < stepnum
							{												
								/* bottom plane 0 */
								theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
								disp_t0.RowAlias(stepnum-l, force1);
								theta.Multx(force1, force2);
								force0a.AddScaled(timestep, force2);	
							}
						}
						else	// normalized time greater than critical value
						{	
							for (int l = 0; l < fNumstep_crit; l++)
							{												
								/* bottom plane 0 */
								theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
								disp_t0.RowAlias(fNumstep_crit-l-1, force1);
								theta.Multx(force1, force2);
								force0a.AddScaled(timestep, force2);	
							}
						}
					}
				}
				count++;
			}
		}
		/* add THK force to fTHKDisp - bottom atoms */
		fTHKforce.SetRow(i, force0b);
		
		if (fnumsets == 2)
		{
			/* add THK force to fTHKDisp - top atoms */
			fTHKforce.SetRow(i+fBottomatoms.Length(), force0a);
		}
	}
	//cout << "THKforce = " << fTHKforce << endl;
	return fTHKforce;

}

/* describe the parameters needed by the interface */
void FEManagerT_THK::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FEManagerT_bridging::DefineParameters(list);

	/* time-history kernel parameters */
	list.AddParameter(ParameterT::Integer, "N_crit");
	list.AddParameter(ParameterT::Double, "lattice_parameter");
	list.AddParameter(ParameterT::Word, "theta_file");
	
	// the ghost mapping for the coupling matrix specified in a file
	list.AddParameter(ParameterT::Word, "ghostmap_file");
}

/* information about subordinate parameter lists */
void FEManagerT_THK::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	FEManagerT_bridging::DefineSubs(sub_list);

	/* nodes affected by THK boundary conditions */
	sub_list.AddSub("THK_nodes_ID_list");
	
}

/* accept parameter list */
void FEManagerT_THK::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	FEManagerT_bridging::TakeParameterList(list);
	
	/* extract THK parameters */
	fNcrit = list.GetParameter("N_crit");
	fLatticeParameter = list.GetParameter("lattice_parameter");
	fThetaFile = list.GetParameter("theta_file");
	fThetaFile.ToNativePathName();
	StringT path;
	path.FilePath(InputFile());
	fThetaFile.Prepend(path);
	
	// extract matrix for use when doing ghostoffmap in BridgingScaleManagerT
	fGhostMapFile = list.GetParameter("ghostmap_file");
	fGhostMapFile.ToNativePathName();
	StringT filepath;
	filepath.FilePath(InputFile());
	fGhostMapFile.Prepend(path);

	/* nodes affected by THK boundary conditions */
	const ParameterListT& id_list = list.GetList("THK_nodes_ID_list");
	StringListT::Extract(id_list, fTHKNodes);

	/* initialize */
	int nsd = NodeManager()->NumSD();
	if (nsd == 2)
		Initialize2D();
	else if (nsd == 3)
		Initialize3D();
	else
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize", 
			"%d dimensions not unsupported", nsd);
}

/*************************************************************************
 * Private
 *************************************************************************/

/* compute theta tables for 2D disp/force formulation */
void FEManagerT_THK::ComputeThetaTables2D(const StringT& data_file)
{
	const char caller[] = "FEManagerT_THK::ComputeThetaTables2D";
	ifstreamT data(data_file);
	if (!data.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
		
	/* dimensions */
	double pi = acos(-1.0);
	int n_neighbor = fNcrit + 1;  // number of neighbors to calculate thetas for
	int n_sum, nsteps;       // number of fourier coefficients to use in calculation of theta
	data >> n_sum;
	
	/* dimension work space */
	const TimeManagerT* time_manager = TimeManager();
	fN_times = time_manager->NumberOfSteps();
	double tstep = time_manager->TimeStep(); // timestep
	double totaltime = fN_times * tstep; // total time of simulation
	
	/* determine correct loop time for theta and time history variables */
	double looptime = 0.0;
	if (totaltime <= 3.0)  // assume that t_crit = 3.0 in normalized time
		looptime = totaltime;
	else
		looptime = 3.0;   // need to normalize this time
	fNumstep_crit = int(3.0/tstep) + 1;	// currently assuming normalized t_crit = 3 
										// need to actually compute using LJ parameters
	/* determine correct number of timesteps to store for theta and history variables */
	if (fN_times+1 <= fNumstep_crit)	// WATCH +1 FACTOR!!!
		nsteps = fN_times+1;  // watch +1
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

/* compute theta tables for 3D disp/disp formulation */
void FEManagerT_THK::ComputeThetaTables3D(const StringT& data_file)
{
	const char caller[] = "FEManagerT_THK::ComputeThetaTables3D";
	ifstreamT data(data_file);
	if (!data.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
						
	/* dimensions */
	double pi = acos(-1.0);
	int n_neighbor = (2 * fNcrit + 1) * (2 * fNcrit + 1);  // number of neighbors to calculate thetas for
	int n_sum, nsteps;       // number of fourier coefficients to use in calculation of theta
	data >> n_sum;
	
	/* dimension work space */
	const TimeManagerT* time_manager = TimeManager();
	fN_times = time_manager->NumberOfSteps();
	double tstep = time_manager->TimeStep(); // timestep
	double totaltime = fN_times * tstep; // total time of simulation

	/* determine correct loop time for theta and time history variables */
	double looptime = 0.0;
	if (totaltime <= 2.0)  // assume that t_crit = 12.0 in normalized time
		looptime = totaltime;
	else
		looptime = 2.0;   // need to normalize this time

	fNumstep_crit = int((2.0/tstep) + 0.5) + 1;	// currently assuming normalized t_crit = 12.0 
	// DEF note: the 0.5 was added on to ensure that the int chop gets the correct number (to avoid roundoff issues)
	
	/* determine correct number of timesteps to store for theta and history variables */
	if (fN_times+1 <= fNumstep_crit)	// WATCH +1 FACTOR!!!
		nsteps = fN_times+1;  // watch +1
	else
		nsteps = fNumstep_crit;

	iArrayT crit(fNumstep_crit);
	fShift.Dimension(fNumstep_crit-1);
	crit.SetValueToPosition();
	fShift.CopyPart(0, crit, 1, fNumstep_crit-1);	// Used to shift displacement history when time > t_crit - MAY NEED TO FIX!
	
	fThetaTableB.Dimension(n_neighbor); //(n_neighbor)
	
	/* dimension boundary atom history */
	fHistoryTableb.Dimension(fBottomatoms.Length());
	
	
	if (fnumsets == 2) 
	{ 	
		fThetaTableT.Dimension(n_neighbor);	//(n_neighbor)
		
		// dimension boundary atom history
		fHistoryTablet.Dimension(fTopatoms.Length());
		
		for (int i = 0; i < fBottomatoms.Length(); i++)
		{
			fHistoryTablet[i].Dimension(nsteps, 3);  // used to be fN_times+1,2
		}
		
		// dimension theta table data structure
		for (int i = 0; i < fThetaTableB.Length(); i++)
		{
			fThetaTableT[i].Dimension(nsteps, 3*3);  // used to be fN_times+1,2*2
		}
		
	}
	
	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		fHistoryTableb[i].Dimension(nsteps, 3);	// used to be fN_times+1,2
	}
	
	// dimension theta table data structure */
	for (int i = 0; i < fThetaTableB.Length(); i++)
	{
		fThetaTableB[i].Dimension(nsteps, 3*3);  // used to be fN_times+1,2*2
	}
	
	/* read in data tables for multiple thetas */
	dArrayT row;
	ArrayT<dArray2DT> data_table(n_neighbor);
	for (int i = 0; i < data_table.Length(); i++)
	{
		/* dimension */
		dArray2DT& n_table = data_table[i];
		n_table.Dimension(n_sum, 3*3);
				
		/* read */
		for (int j = 0; j < n_table.MajorDim(); j++)
		{
			int junk;
			data >> junk >> junk >> junk;
			
			/* each row is b^T */
			n_table.RowAlias(j, row);
			data >> row;
		}
	}
		
	dMatrixT theta1(3), temptheta(3), transform(3);
	transform.Identity(1.0);
	transform(2,2) = -1.0; // Theta(bottom)=T*Theta(top)*T, T=transform

	/* compute theta's - top/bottom planes related by transform matrix above */
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
				temptheta.AddScaled(sin((k+1)*pi*j/2.2), theta_i(k)); // tmax = 2.2 now
			}

			/* add temptheta into fThetaTable */
			if (fnumsets == 2) // top theta table -> may wish to pull out to save time
			{
				fThetaTableT[i].SetRow(count, temptheta);
			}
		
			/* compute bottom thetatables via transform matrix */
			theta1.MultABC(transform, temptheta, transform);
			fThetaTableB[i].SetRow(count, theta1);
			count++;
		}
	}
}

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
