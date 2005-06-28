/* $Id: FEManagerT_THK.cpp,v 1.23 2005-06-28 14:45:39 d-farrell2 Exp $ */

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
#include "ParameterContainerT.h"

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
	int nsd = model->NumDimensions();
	
	// read other parameters and initialize data
	fNeighbors = 2 * fNcrit + 1;   // maximum number of neighbors per atom in 2D

	/* read node set indexes */
	fnumsets = fTHKNodes.Length();		// number of MD THK boundary node sets
 
	/* collect sets - currently assuming exactly 2 MD THK boundaries in 2D */
	/* further assumes bottom is first node set, top is second node set */
#pragma message("FEManagerT_THK::Initialize2D, only set up for consistent input with 3D. -- talk to dave")
	
	// collect sets, set up some arrays (these have repeats in general)
	if (fnumsets >= 1)
	{
		if (fnumsets > 4)
		{
#pragma message("The formulation used here is not good if you have a square lattice like structure as a boundary, needs to be generalized further")
			cout << "THe formulation used here is only good when in neighboring layers, there is not an atom in the same location" << endl;
		}
		fbound_set_atoms.Dimension(fnumsets);
		fbound_neighbor_atoms.Dimension(fnumsets);
		for (int i = 0; i < fnumsets; i++)
		{
			// collect sets - minimum 1
			// put each set into array of integer arrays
			fbound_set_atoms[i] = model->NodeSet(fTHKNodes[i]);
			
			// dimension array of 2D arrays to hold in-plane neighbors ( host and other sets as well)
			// set the value to -1, for no neighbor
			fbound_neighbor_atoms[i].Dimension(fbound_set_atoms[i].Length(), fNeighbors);
			fbound_neighbor_atoms[i] = -1;
		}
	}
	else
	{
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D - # of THK BC sets less than 1"); // this is the only check on fnumsets
	}
	
	// figure out how many boundary atoms there are (since this is based on sets with repeats, this is corrected after the special atom info is read in)
	ftotal_b_atoms = 0;
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		ftotal_b_atoms += fbound_set_atoms[isetcount].Length(); 
	}
	
	// Parse the special atom information file
	ifstreamT data(fSpecAtomFile);
	if (!data.is_open())
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D", "file not found: %s", fSpecAtomFile.Pointer());
	
	// now grab data from file
	
	// first entry in file is number of special atoms
	data >> fnum_spec_atoms;
	
	if (fnum_spec_atoms > 0) // if there are special atoms
	{
		// now read in the data, need to set up some stuff first
		fSpecAtomInfo.Dimension(fnum_spec_atoms);
		fSpecAtomID.Dimension(fnum_spec_atoms);
		int tempatomnum, tempnumsets, temphostset;
		iArrayT tempspecinfo;
		
		int tempspeccounter = 0;
		
		// grab the special atom data
		for (int i = 0; i < fnum_spec_atoms; i++ ) // for all of the special atoms
		{
			// read in atom #, # of boundaries it is a member of (2 -> fNumSets), host boundary # (1 -> fNumSets) (set the atom is found in in fbound_set_atoms)
			data >> tempatomnum >> tempnumsets >> temphostset ;
			if (tempnumsets < 2)
				ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D", "Special Atom %d, number of member sets should be >= 2, not %d", tempatomnum, tempnumsets);
			
			tempspecinfo.Dimension(tempnumsets+2); // dimension the temp array to hold the information
			tempspecinfo = -1;	// reset the array to -1 for all entries
			tempspecinfo[0] = tempatomnum;
			tempspecinfo[1] = tempnumsets;
			tempspecinfo[2] = temphostset;
			
			// now read in remaining sets
			for (int j = 0; j < (tempnumsets-1); j++)
			{
				data >> tempspecinfo[j+3]; // j+3 to make sure alignment is correct
			}
			
			// now place temporary array in the permanant array
			fSpecAtomInfo[i] = tempspecinfo;		// this contains ALL of the information needed to deal with special atoms, but may not be needed in its entirety
			fSpecAtomID[i] = tempatomnum-1;			// this offsets the value to the internal atom numbering 0-(NN-1) NN = # nodes
			tempspeccounter +=  (tempnumsets-1);	// counter to keep record of how many extra entries there are (tempnumsets-1 is since there should be 1 entry for each)
		}
		
		// now correct ftotal_b_atoms to remove multiple counting of special atoms
		ftotal_b_atoms -= tempspeccounter;

	}
	
	// dimension array for THK force (sized for 1 entry per atom on boundary, no repeats)
	fTHKforce.Dimension(ftotal_b_atoms, 2);
#pragma message("Dave set up the initialize2D to be basically the same as 3D. this is for consistency, expense not big issue since only done once.")

	// do neighbor search - only done once as of original implmententation -- by assumptions of method, only needs to be done once
	// this looks for neighboring atoms in the boundary plane, which are to be considered in the summation mentioned in eqn 41 in hspark's 2D paper
	NodeManagerT* node = FEManagerT::NodeManager();				// get the node manager
	const dArray2DT& initcoords = node->InitialCoordinates();	// get initial coords for all atoms
	iArrayT temp_atom(initcoords.MajorDim());					// set up another temp integer array with length = # atoms
	temp_atom.SetValueToPosition();								// initialize the above array to a local numbering scheme
	
	/* configure search grid - CURRENTLY SEARCHING ONLY NON-IMAGE ATOMS (REAL+GHOST) */
	iGridManagerT grid(10, 100, initcoords, &temp_atom);		// set up search grid based on initial coordinates of real + ghost atoms
	grid.Reset();
	dArrayT acoord1, acoord2, ncoord1, ncoord2;					// declare some more temp arrays (maybe we can kill some of these off, save memory?)
	double checkdot, checkrad, mag1, xmag, ymag;			// declare some temp variables
	int counter;
	dArrayT facenormal;											// nodeset normal (error = 0)											
	facenormal = 0.0;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{ 
		
		// get the normal of plane of boundary atoms (provided in input file)
		facenormal = fTHK_normals[isetcount];
		if (facenormal == 0.0)
			ExceptionT::GeneralFail("FEManagerT_THK::Initialize2D", "THK Boundary face normal not defined");
		
		// find in-plane neighbor atoms, based on results of search (2 level search, coarse and fine)
		// first do regular atoms and 'special' atoms with respect to each set to which it belongs
		for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
		{
			initcoords.RowAlias(fbound_set_atoms[isetcount][i], acoord1);		// alias the coordinates of the ith real boundary atom to acoord1
			
			// access the candidate neighbor points
			const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord1.Pointer(), (fNcrit+1)*1.05*sqrt(2.0)*fLatticeParameter); // will work for FCC, BCC, SC 
																																// DEF added (fNCrit+1) to scale the search area
			for (int j = 0; j < hitsa.Length(); j++)		// loop over the candidate neighbors
			{			
				// distance between atom and candidate neighbor
				ncoord1.Alias(nsd, hitsa[j].Coords());		// alias the coordinates of the candidate neighbor to ncoord1 (perhaps change syntax?)
				xmag = ncoord1[0]-acoord1[0];				// get the distance between the atom and candidate neighbor point in each direction
				ymag = ncoord1[1]-acoord1[1];
				mag1 = sqrt(xmag*xmag+ymag*ymag);	// calculate the point to point distance between the neigbor and atom
				
				// now do neighbor search first determine which of the candidates is in the same plane (look to see if facenormal & normalized vector between atoms is ~ +/- tol)
				// once the atoms are determined to be in the same plane, look at the distance, see if it is within the desired range.
				// this will work for arbitrary plane orientations - but lattice parameter is the in-plane parameter
							
				// Get in-plane neighbors
				if (mag1 > 1.0e-8) // if mag1 is close enough to be zero -> this would blow up -> but should be zero (allows self checking)
					checkdot = (1/mag1) *((facenormal[0] * xmag) + (facenormal[1] * ymag)); // compute dot product of normal and candidate vector
				else
					checkdot = 0.0;
				
				if (fabs(checkdot) < tol )	// if the candidate neighbor is within the plane, check it
				{																					 					
					// do fine check, store matches in array for THK BC application
					
					// fNCrit denotes neighborshell in the boundary plane
					// if fNCrit = 0, self is only one considered in BC summation
					counter = 0;
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						// check distance against position in neighbor shell
						checkrad = sqrt((fLatticeParameter*l)*(fLatticeParameter*l));
						if ( fabs(mag1-checkrad) < tol)
						{																				
							fbound_neighbor_atoms[isetcount](i,counter) = hitsa[j].Tag();
							counter++;	// counter increments to make sure atoms are in right order									
						}
					}																				
				}																					
			}																					
		}
	}
	
	// obtain the ghost atom properties properties map (to turn off the atoms when computing the FEM force
	DoGhostMap();
	
	
	/* compute theta tables */
	ComputeThetaTables();
}

/* Bridging Scale 3D Initialization */
void FEManagerT_THK::Initialize3D(void)
{
	/* Implement 3D version of initialize here */
	/* Find neighbors of top 2 planes of atoms */
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();

	// get information on THK BC sets

	// set number of neighbor atoms to consider for the THK BC (including the boundary atom)
	int num_neighborsx = 2 * fNcrit + 1;	// assume same number of neighbors in both x and y directions 
	int num_neighborsy = 2 * fNcrit + 1;
	fNeighbors = num_neighborsx * num_neighborsy; // total number of neighbors to consider in the sum on the THK terms
												  // see eqn 41 in hspark's 2D paper (for now, assume same for all sets)
	// read node set indices
	fnumsets = fTHKNodes.Length();		// number of MD THK boundary node sets
	
	// collect sets, set up some arrays (these have repeats)
	if (fnumsets >= 1)
	{
		if (fnumsets > 6)
		{
#pragma message("The formulation used here is not good if you have a simple cubic like structure as a boundary, needs to be generalized further")
			cout << "THe formulation used here is only good when in neighboring layers, there is not an atom in the same l,m location" << endl;
		}
		fbound_set_atoms.Dimension(fnumsets);
		fbound_neighbor_atoms.Dimension(fnumsets);
		for (int i = 0; i < fnumsets; i++)
		{
			// collect sets - no maximum, minimum 1
			// put each set into array of integer arrays
			fbound_set_atoms[i] = model->NodeSet(fTHKNodes[i]);
			
			// dimension array of 2D arrays to hold in-plane neighbors ( host and other sets as well)
			// set the value to -1, for no neighbor
			fbound_neighbor_atoms[i].Dimension(fbound_set_atoms[i].Length(), fNeighbors);
			fbound_neighbor_atoms[i] = -1;
		}
	}
	else
	{
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D - # of THK BC sets less than 1, check inputs"); // this is the only check on fnumsets
	}
	
	// figure out how many boundary atoms there are (since this is based on sets with repeats, this is corrected after the special atom info is read in)
	ftotal_b_atoms = 0;
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		ftotal_b_atoms += fbound_set_atoms[isetcount].Length(); 
	}
	
///// May wish to move this to another method
	
	// Parse the special atom information file
	ifstreamT data(fSpecAtomFile);
	if (!data.is_open())
		ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D", "file not found: %s", fSpecAtomFile.Pointer());
	
	// now grab data from file
	
	// first entry in file is number of special atoms
	data >> fnum_spec_atoms;
	
	if (fnum_spec_atoms > 0) // if there are special atoms
	{
		// now read in the data, need to set up some stuff first
		fSpecAtomInfo.Dimension(fnum_spec_atoms);
		fSpecAtomID.Dimension(fnum_spec_atoms);
		int tempatomnum, tempnumsets, temphostset;
		iArrayT tempspecinfo;
		
		int tempspeccounter = 0;
		
		// grab the special atom data
		for (int i = 0; i < fnum_spec_atoms; i++ ) // for all of the special atoms
		{
			// read in atom #, # of boundaries it is a member of (2 -> fNumSets), host boundary # (1 -> fNumSets) (set the atom is found in in fbound_set_atoms)
			data >> tempatomnum >> tempnumsets >> temphostset ;
			if (tempnumsets < 2)
				ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D", "Special Atom %d, number of member sets should be >= 2, not %d", tempatomnum, tempnumsets);
			
			tempspecinfo.Dimension(tempnumsets+2); // dimension the temp array to hold the information
			tempspecinfo = -1;	// reset the array to -1 for all entries
			tempspecinfo[0] = tempatomnum;
			tempspecinfo[1] = tempnumsets;
			tempspecinfo[2] = temphostset;
			
			// now read in remaining sets
			for (int j = 0; j < (tempnumsets-1); j++)
			{
				data >> tempspecinfo[j+3]; // j+3 to make sure alignment is correct
			}
			
			// now place temporary array in the permanant array
			fSpecAtomInfo[i] = tempspecinfo;		// this contains ALL of the information needed to deal with special atoms, but may not be needed in its entirety
			fSpecAtomID[i] = tempatomnum-1;			// this offsets the value to the internal atom numbering 0-(NN-1) NN = # nodes
			tempspeccounter +=  (tempnumsets-1);	// counter to keep record of how many extra entries there are (tempnumsets-1 is since there should be 1 entry for each)
		}
		
		// now correct ftotal_b_atoms to remove multiple counting of special atoms
		ftotal_b_atoms -= tempspeccounter;

	}
/////
		
	// dimension array for THK force (sized for 1 entry per atom on boundary, no repeats)
	fTHKforce.Dimension(ftotal_b_atoms, 3);


	// do neighbor search - only done once as of original implmententation -- by assumptions of method, only needs to be done once
	// this looks for neighboring atoms in the boundary plane, which are to be considered in the summation mentioned in eqn 41 in hspark's 2D paper
	NodeManagerT* node = FEManagerT::NodeManager();				// get the node manager
	const dArray2DT& initcoords = node->InitialCoordinates();	// get initial coords for all atoms
	iArrayT temp_atom(initcoords.MajorDim());					// set up another temp integer array with length = # atoms
	temp_atom.SetValueToPosition();								// initialize the above array to a local numbering scheme
	
	/* configure search grid - CURRENTLY SEARCHING ONLY NON-IMAGE ATOMS (REAL+GHOST) */
	iGridManagerT grid(10, 100, initcoords, &temp_atom);		// set up search grid based on initial coordinates of real + ghost atoms
	grid.Reset();
	dArrayT acoord1, acoord2, ncoord1, ncoord2;					// declare some more temp arrays (maybe we can kill some of these off, save memory?)
	double checkdot, checkrad, mag1, xmag, ymag, zmag;			// declare some temp variables
	int counter;
	dArrayT facenormal;											// nodeset normal (error = 0)											
	facenormal = 0.0;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{ 
		
		// get the normal of plane of boundary atoms (provided in input file)
		facenormal = fTHK_normals[isetcount];
		if (facenormal == 0.0)
			ExceptionT::GeneralFail("FEManagerT_THK::Initialize3D", "THK Boundary face normal not defined");
		
		// find in-plane neighbor atoms, based on results of search (2 level search, coarse and fine)
		// first do regular atoms and 'special' atoms with respect to each set to which it belongs
		for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
		{
			initcoords.RowAlias(fbound_set_atoms[isetcount][i], acoord1);		// alias the coordinates of the ith real boundary atom to acoord1
			
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
				// this will work for arbitrary plane orientations - but lattice parameter is the in-plane parameter
							
				// Get in-plane neighbors
				if (mag1 > 1.0e-8) // if mag1 is close enough to be zero -> this would blow up -> but should be zero (allows self checking)
					checkdot = (1/mag1) *((facenormal[0] * xmag) + (facenormal[1] * ymag) + (facenormal[2] * zmag)); // compute dot product of normal and candidate vector
				else
					checkdot = 0.0;
				
				if (fabs(checkdot) < tol )	// if the candidate neighbor is within the plane, check it
				{																					 					
					// do fine check, store matches in array for THK BC application
					
					// fNCrit denotes neighborshell in the boundary plane
					// if fNCrit = 0, self is only one considered in BC summation
					counter = 0;
					for (int k = -fNcrit; k <= fNcrit; k++)
					{
						for (int l = -fNcrit; l<= fNcrit; l++)
						{
							// check distance against position in neighbor shell
							checkrad = sqrt((fLatticeParameter*k)*(fLatticeParameter*k) + (fLatticeParameter*l)*(fLatticeParameter*l));
							if ( fabs(mag1-checkrad) < tol)
							{																				
								fbound_neighbor_atoms[isetcount](i,counter) = hitsa[j].Tag();
								counter++;		// counter increments to make sure atoms are in right order									
							}
						}																			
					}																				
				}																					
			}																					
		}
	}
	
	// obtain the ghost atom properties properties map (to turn off the atoms when computing the FEM force
	DoGhostMap();
	
	// have boundary atom THK neighbor list, have properties mapping for calculation of forces
	// now Compute Theta tables
	ComputeThetaTables();
}


// 2D/3D MD/THK Initialization
void FEManagerT_THK::InitializeMDTHK(void)
{
#pragma message("The formulation used here is not good if you have a square lattice like structure as a boundary (for non-nearest neighbors), needs to be generalized further")	
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
	
	// read other parameters and initialize data
	if (nsd == 2)
		fNeighbors = 2 * fNcrit + 1;	// maximum number of neighbors per atom in 2D
	else if (nsd == 3)
		fNeighbors = (2 * fNcrit + 1)*(2 * fNcrit + 1);	// maximum number of neighbors per atom in 3D
	else
		ExceptionT::GeneralFail("FEManagerT_THK::InitializeMDTHK", "%d dimensions not implemented", nsd);

	/* read node set indexes */
	fnumsets = fTHKNodes.Length();				// number of MD THK boundary node sets
	int numsets_temp = fTHKGhostNodes.Length();	// number of MD THK ghost node sets
	
	if (fnumsets != numsets_temp)
		ExceptionT::GeneralFail("FEManagerT_THK::InitializeMDTHK - # of THK BC sets does not match # of THK Ghost atom Sets");
	
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
		ExceptionT::GeneralFail("FEManagerT_THK::InitializeMDTHK - # of THK BC sets less than 1"); // this is the only check on fnumsets
	}
	
	// figure out how many ghost atoms there are (probably could just ask the nodemanager)
	ftotal_b_atoms = 0;
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		ftotal_b_atoms += fghost_set_atoms[isetcount].Length(); 
	}
	
	// dimension array for THK force (sized for 1 entry per atom on boundary, no repeats)
	fTHKdisp.Dimension(ftotal_b_atoms, nsd);

	// do neighbor search - only done once as of original implmententation -- by assumptions of method, only needs to be done once
	// this looks for neighboring atoms in the boundary plane (_not_ the ghost plane)
	if (nsd == 2)
		DoNeighSearch2D();
	else if (nsd == 3)
		DoNeighSearch3D();
	else
		ExceptionT::GeneralFail("FEManagerT_THK::InitializeMDTHK", "%d dimensions not implemented", nsd);
	
	// compute theta tables
	ComputeThetaTables();
}

/* return iArrayT of boundary and ghost atom numbers - 3D version */
const iArrayT& FEManagerT_THK::InterpolationNodes3D(void)
{
	const iArrayT& ghost = GhostNodes();
	fInterpolationNodes.Dimension(ghost.Length()+ftotal_b_atoms);
	fInterpolationNodes.CopyIn(0,ghost);	// copy in ghost atoms first
	
	
	// loop over each boundary, put boundary atoms into the array(make it so it has no repeats)
	
	if (fnum_spec_atoms > 0)	// if there are special atoms, have to monkey with things a bit
	{
		int offset = -1;
		int specpos;
		for (int isetcount = 0; isetcount < fnumsets; isetcount++)
		{
			for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it
					{
						offset++;	// increment the counter
						fInterpolationNodes[(ghost.Length() + offset)] = fbound_set_atoms[isetcount][i];	// place the value
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set, dont count it
					{
						continue;	// go to the next value, 
					}
				}
				else	// not special
				{
					offset++;	// increment the counter
					fInterpolationNodes[(ghost.Length() + offset)] = fbound_set_atoms[isetcount][i];	// place the value
				}
			}
			
		}		
	}
	else	// no changes to original implementation
	{
		int offset = 0;
		for (int isetcount = 0; isetcount < fnumsets; isetcount++)
		{
			fInterpolationNodes.CopyIn((ghost.Length() + offset), fbound_set_atoms[isetcount]);	// then copy in set atoms
			offset += fbound_set_atoms[isetcount].Length();
		}		
	}

	
	return fInterpolationNodes;	
}

/* return iArrayT of boundary and ghost atom numbers - 2D version */
const iArrayT& FEManagerT_THK::InterpolationNodes2D(void)
{
	const iArrayT& ghost = GhostNodes();
	fInterpolationNodes.Dimension(ghost.Length()+ftotal_b_atoms);
	fInterpolationNodes.CopyIn(0,ghost);	// copy in ghost atoms first
	
	
	// loop over each boundary, put boundary atoms into the array(make it so it has no repeats)
	
	if (fnum_spec_atoms > 0)	// if there are special atoms, have to monkey with things a bit
	{
		int offset = -1;
		int specpos;
		for (int isetcount = 0; isetcount < fnumsets; isetcount++)
		{
			for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it
					{
						offset++;	// increment the counter
						fInterpolationNodes[(ghost.Length() + offset)] = fbound_set_atoms[isetcount][i];	// place the value
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set, dont count it
					{
						continue;	// go to the next value, 
					}
				}
				else	// not special
				{
					offset++;	// increment the counter
					fInterpolationNodes[(ghost.Length() + offset)] = fbound_set_atoms[isetcount][i];	// place the value
				}
			}
			
		}		
	}
	else	// no changes to original implementation
	{
		int offset = 0;
		for (int isetcount = 0; isetcount < fnumsets; isetcount++)
		{
			fInterpolationNodes.CopyIn((ghost.Length() + offset), fbound_set_atoms[isetcount]);	// then copy in set atoms
			offset += fbound_set_atoms[isetcount].Length();
		}		
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

/* calculate external force on MD boundary atoms for 2D disp/force formulation */
const dArray2DT& FEManagerT_THK::THKForce2D(const StringT& bridging_field, const dArray2DT& badisp)
{
	// badisp is in format with each set sequentially in vector form
	fTHKforce = 0.0;  

	// dimension some local arrays (these get reused and re defined constantly here)
	dArrayT atomdisp, femdisp(2), diff(2);
	
	/* access the actual MD displacements */
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);
	dArray2DT mddisp = (*atomfield)[0];

	const int stepnum = FEManagerT::StepNumber();  // to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();  // delta t_md

	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, 2);	// copy all rows of history except last row

	// Calculate q - ubar for each set
	
	int dispcounter = -1; // keeps track of place in badisp array
	int specpos;
	iArrayT thk_force_pos_spec(fnum_spec_atoms);
	thk_force_pos_spec = -1;
	int pos_temp = -1;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// now go through the set
		for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
		{
			if (fnum_spec_atoms > 0)	// if there are special atoms, have to use some special considerations (perhaps this can be done better?)
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it ( zero everything, since it is first time seen)
					{
						dispcounter++;	// increment the counter
						pos_temp = dispcounter;
						thk_force_pos_spec[specpos] = dispcounter;	// note the position
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set, access previous record
					{
						pos_temp = thk_force_pos_spec[specpos];					
					}
				}
				else	// not special
				{
					dispcounter++;	// increment the counter
					pos_temp = dispcounter;
				}
			}
			else	// no change
			{
				dispcounter++;	// increment the counter
				pos_temp = dispcounter;	// increment the counter
			}			
			
			/* non-shift case */
			if (stepnum < fNumstep_crit)   
			{
				// for the boundary plane
				mddisp.RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
				badisp.RowAlias(pos_temp, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				fHistoryTable[isetcount][i].SetRow(stepnum, diff);
			}
			else	// t > t_crit
			{
				// for the boundary plane
				mddisp.RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
				badisp.RowAlias(pos_temp, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				shift.RowCollect(fShift, fHistoryTable[isetcount][i]);
				fHistoryTable[isetcount][i].BlockRowCopyAt(shift, 0);
				fHistoryTable[isetcount][i].SetRow(fNumstep_crit-1, diff);
			}
		}
	}

	dMatrixT theta(2);
	dArrayT force1(2), force2(2), force0a(2), force0b(2);
	InverseMapT setnodes;
	int counter2 = -1;
	specpos = -1;
	thk_force_pos_spec = -1;

	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// set inverse maps  - what does this do?? 
		setnodes.SetMap(fbound_set_atoms[isetcount]);	// Set global to local map for bottom atoms plane 0
		int count, dex;

		/* calculate THK force for the boundary plane of atoms */
		for (int i = 0; i < fbound_neighbor_atoms[isetcount].MajorDim(); i++)
		{
			if (fnum_spec_atoms > 0)	// if there are special atoms, have to use some special considerations
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it ( zero everything, since it is first time seen)
					{
						counter2++;	// increment the counter
						thk_force_pos_spec[specpos] = counter2;	// note the position
						force0a = 0.0;	// zero the forces
						force0b = 0.0;
						
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set
					{
						force0a = fTHKforce[thk_force_pos_spec[specpos]];	// record the value of the force already present
					}
				}
				else	// not special
				{
					counter2++;	// increment the counter
					force0a = 0.0;	// zero the forces
					force0b = 0.0;
				}
			}
			else	// no change
			{
				counter2++;	// increment the counter
				force0a = 0.0;	// zero the forces
				force0b = 0.0;
			}
			
			count = 0;
			for (int j = 0; j < 2*fNcrit+1; j++)
			{
				if (fbound_neighbor_atoms[isetcount](i,count) == -1)	// if no neighbor
					int temp = 0;
				else
				{
					/* bottom plane 0 */
					dex = setnodes.Map(fbound_neighbor_atoms[isetcount](i,count));
					const dArray2DT& disp_temp = fHistoryTable[isetcount][dex];
					const dArray2DT& theta_temp = fThetaTable_array[isetcount][fNeighbors-1-count]; 	// start count at 8, go backwards 
					
					/* calculate fine scale THK disp using theta here */
					if (stepnum < fNumstep_crit)  // < fNumstep_crit 
					{
						for (int l = 0; l < stepnum; l++)
						{												
							theta.Alias(2,2,theta_temp(l));
							disp_temp.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{												
							theta.Alias(2,2,theta_temp(l));
							disp_temp.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	
						}
					}
				}
				count++;
			}
			// THK Force
			if (fnum_spec_atoms > 0)	// if there are special atoms, have to use some special considerations
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it ( zero everything, since it is first time seen)
					{
						fTHKforce.SetRow(counter2, force0b);
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set
					{
						force0b += force0a;	// simply add the just calculated force to the previous (this may need to change, since simple sum may not be exact)
						fTHKforce.SetRow(thk_force_pos_spec[specpos], force0b);	// then set the value
					}
				}
				else	// not special
				{
					fTHKforce.SetRow(counter2, force0b);
				}
			}
			else	// no change
			{
				fTHKforce.SetRow(counter2, force0b);
			}
		}
	}
	//cout << "THKforce = " << fTHKforce << endl;
	return fTHKforce;
}

/* calculate impedance force using 3D disp/force formulation */
const dArray2DT& FEManagerT_THK::THKForce3D(const StringT& bridging_field, const dArray2DT& badisp)
{
	// badisp is in format with each set sequentially in vector form
	fTHKforce = 0.0;

	// dimension some local arrays (these get reused and re defined constantly here)
	dArrayT atomdisp, femdisp(3), diff(3);
	
	/* access the actual MD displacements */
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);
	dArray2DT mddisp = (*atomfield)[0];

	const int stepnum = FEManagerT::StepNumber();  // to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();  // delta t_md

	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, 3);	// copy all rows of history except last row 

	// Calculate q - ubar for each set
	
	int dispcounter = -1; // keeps track of place in badisp array
	int specpos;
	iArrayT thk_force_pos_spec(fnum_spec_atoms);
	thk_force_pos_spec = -1;
	int pos_temp = -1;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// now go through the set
		for (int i = 0; i < fbound_set_atoms[isetcount].Length(); i++)
		{
			if (fnum_spec_atoms > 0)	// if there are special atoms, have to use some special considerations (perhaps this can be done better?)
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it ( zero everything, since it is first time seen)
					{
						dispcounter++;	// increment the counter
						pos_temp = dispcounter;
						thk_force_pos_spec[specpos] = dispcounter;	// note the position
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set, access previous record
					{
						pos_temp = thk_force_pos_spec[specpos];					
					}
				}
				else	// not special
				{
					dispcounter++;	// increment the counter
					pos_temp = dispcounter;
				}
			}
			else	// no change
			{
				dispcounter++;	// increment the counter
				pos_temp = dispcounter;	// increment the counter
			}			
			
			/* non-shift case */
			if (stepnum < fNumstep_crit)   
			{
				// for the boundary plane
				mddisp.RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
				badisp.RowAlias(pos_temp, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				fHistoryTable[isetcount][i].SetRow(stepnum, diff);
			}
			else	// t > t_crit
			{
				// for the boundary plane
				mddisp.RowAlias(fbound_set_atoms[isetcount][i], atomdisp);
				badisp.RowAlias(pos_temp, femdisp);
				diff.DiffOf(atomdisp, femdisp);
				shift.RowCollect(fShift, fHistoryTable[isetcount][i]);
				fHistoryTable[isetcount][i].BlockRowCopyAt(shift, 0);
				fHistoryTable[isetcount][i].SetRow(fNumstep_crit-1, diff);
			}
		}
	}
	
	dMatrixT theta(3);
	dArrayT force1(3), force2(3), force0a(3), force0b(3);
	InverseMapT setnodes;
	int counter2 = -1;
	specpos = -1;
	thk_force_pos_spec = -1;
	
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// set inverse maps  - what does this do?? 
		setnodes.SetMap(fbound_set_atoms[isetcount]);	// Set global to local map for bottom atoms plane 0
		int count, dex;
		//dex = setnodes.Map(839);

		/* calculate THK force for the boundary plane of atoms */
		for (int i = 0; i < fbound_neighbor_atoms[isetcount].MajorDim(); i++)
		{
			if (fnum_spec_atoms > 0)	// if there are special atoms, have to use some special considerations
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it ( zero everything, since it is first time seen)
					{
						counter2++;	// increment the counter
						thk_force_pos_spec[specpos] = counter2;	// note the position
						force0a = 0.0;	// zero the forces
						force0b = 0.0;
						
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set
					{
						force0a = fTHKforce[thk_force_pos_spec[specpos]];	// record the value of the force already present
					}
				}
				else	// not special
				{
					counter2++;	// increment the counter
					force0a = 0.0;	// zero the forces
					force0b = 0.0;
				}
			}
			else	// no change
			{
				counter2++;	// increment the counter
				force0a = 0.0;	// zero the forces
				force0b = 0.0;
			}
			
			count = 0;
			for (int j = 0; j < 2*fNcrit+1; j++)
			{
				for (int k = 0; k < 2*fNcrit+1; k++)
				{
					if (fbound_neighbor_atoms[isetcount](i,count) == -1)	// if no neighbor
						int temp = 0;
					else
					{
						/* bottom plane 0 */
						dex = setnodes.Map(fbound_neighbor_atoms[isetcount](i,count));
						const dArray2DT& disp_temp = fHistoryTable[isetcount][dex];
						const dArray2DT& theta_temp = fThetaTable_array[isetcount][fNeighbors-1-count]; 	// start count at 8, go backwards 
						
						/* calculate fine scale THK disp using theta here */
						if (stepnum < fNumstep_crit)  // < fNumstep_crit 
						{
							for (int l = 0; l < stepnum; l++)
							{												
								theta.Alias(3,3,theta_temp(l));
								disp_temp.RowAlias(stepnum-l, force1);
								theta.Multx(force1, force2);
								force0b.AddScaled(timestep, force2);	
							}
						}
						else	// normalized time greater than critical value
						{	
							for (int l = 0; l < fNumstep_crit; l++)
							{												
								theta.Alias(3,3,theta_temp(l));
								disp_temp.RowAlias(fNumstep_crit-l-1, force1);
								theta.Multx(force1, force2);
								force0b.AddScaled(timestep, force2);	
							}
						}
					}
					count++;
				}
			}
			// THK Force
			if (fnum_spec_atoms > 0)	// if there are special atoms, have to use some special considerations
			{
				if (fSpecAtomID.HasValue(fbound_set_atoms[isetcount][i], specpos) == 1)	// if current atom is special
				{
					if ((fSpecAtomInfo[specpos][2] - 1) == isetcount )	// host set, count it ( zero everything, since it is first time seen)
					{
						fTHKforce.SetRow(counter2, force0b);
					}
					else if ((fSpecAtomInfo[specpos][2] - 1) != isetcount) // not host set
					{
						force0b += force0a;	// simply add the just calculated force to the previous (this may need to change, since simple sum may not be exact)
						fTHKforce.SetRow(thk_force_pos_spec[specpos], force0b);	// then set the value
					}
				}
				else	// not special
				{
					fTHKforce.SetRow(counter2, force0b);
				}
			}
			else	// no change
			{
				fTHKforce.SetRow(counter2, force0b);
			}
		}
	}
	//cout << "THKforce = " << fTHKforce << endl;
	return fTHKforce;
}

/*  calculate THK displacement for ghost atoms for 2/3D disp formulation (set up for 3D now)*/
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
	 * 3)	For now, it is set up for a 'fixed' boundary. Once running, put in second THK
	 *		and Linear gradient model choices... perhaps some user inputs/choices
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

	// Calculate q - ubar for each set
	
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

		// calculate the ghost atom displacement u_g
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

/* describe the parameters needed by the interface */
void FEManagerT_THK::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	FEManagerT_bridging::DefineParameters(list);

	/* time-history kernel parameters */
	list.AddParameter(ParameterT::Integer, "N_crit");
	list.AddParameter(ParameterT::Double, "T_crit");	// user can set T-Crit to match theta files
	list.AddParameter(ParameterT::Double, "lattice_parameter");
	list.AddParameter(ParameterT::Double, "interplanar_parameter");	// for finding boundary plane neighbors
	
	ParameterT omega_sys(ParameterT::Double, "Omega_sys");
	omega_sys.SetDefault(1.0);
	list.AddParameter(omega_sys);	// parameter to scale Tcrit and the bn's if the k used in theta file differs
	
	// the ghost mapping for the coupling matrix specified in a file
	list.AddParameter(ParameterT::Word, "ghostmap_file");
	
	// file which contains 'special' atom info (corners, edges)
	list.AddParameter(ParameterT::Word,"THK_special_atom_file");
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
	fTcrit = list.GetParameter("T_crit");
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
	
	// extract the special atom info
	fSpecAtomFile = list.GetParameter("THK_special_atom_file");
	fSpecAtomFile.ToNativePathName();
	fSpecAtomFile.Prepend(path);
	
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
	// found correctly (basically, assumes that _only_ one of the entries in the normal vector is 1)
	
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
	fghostoffmap = 0;
	
	int ghosti = 0; // row number
	int ghostj = 0; // column number
	// now read in the rows and columns of the fghostoffmap 2D array
	for (int ghostk = 0; ghostk < 2*(promap_dim-1); ghostk++)
	{
		data2 >> ghosti >> ghostj ;
		fghostoffmap(ghosti,ghostj) = 1;
	}
}

// compute theta tables for 2D/3D disp/disp or disp/force formulation (doesn't matter, its all the same)
void FEManagerT_THK::ComputeThetaTables(void)
{	
	const char caller[] = "FEManagerT_THK::ComputeThetaTablesMDTHK";
	
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();
		
	/* dimensions */
	double pi = acos(-1.0);
	int n_sum, nsteps;       // number of fourier coefficients to use in calculation of theta
	
	/* dimension work space */
	const TimeManagerT* time_manager = TimeManager();	// get the time manager
	fN_times = time_manager->NumberOfSteps();		 	// get the number of timesteps
	double tstep = time_manager->TimeStep(); 		 	// get the timestep size
	double totaltime = fN_times * tstep; 				// total time of simulation
	
	/* determine correct loop time for theta and time history variables */
	double looptime = 0.0;
	if (totaltime <= .75*(fTcrit/fOmega_sys))  // need to store/calculate up to this time(constant to avoid uglyness at end of period)
		looptime = totaltime;
	else
		looptime = .75*(fTcrit/fOmega_sys);   // need to store/calculate up to this time

	fNumstep_crit = int(.75*(fTcrit/(tstep*fOmega_sys)) + 0.5) + 1;	// determine number of steps to t_crit
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
		
	for (int isetcount = 0; isetcount < fnumsets; isetcount++)
	{
		// Open the file with the fourier coeffs
		const StringT& data_file = fThetaFile_array[isetcount];
		ifstreamT data(data_file);
		if (!data.is_open())
			ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
		
		data >> n_sum;
		
		// Dimension some arrays of arrays to hold information (dimension to # neighbors, # boundary atoms)
		fThetaTable_array[isetcount].Dimension(fNeighbors);
		fHistoryTable[isetcount].Dimension(fbound_set_atoms[isetcount].Length());
		
		// Now dimension all the way down (this has repeats)
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
				int junk;
				if (nsd == 2)
					data >> junk >> junk ; // first 2 cols in file are not needed
				else if (nsd == 3 )
					data >> junk >> junk >> junk ; // first 3 cols in file are not needed
				
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
				
				for (int k = 0; k < n_theta; k++)  // can truncate this summation
				{
					/* Extract each row of coefficients */
					temptheta.AddScaled(fOmega_sys*sin((k+1)*pi*j/(fTcrit/fOmega_sys)), theta_i(k));
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

#endif /* BRIDGING_ELEMENT && BRIDGING_ELEMENT_DEV */
