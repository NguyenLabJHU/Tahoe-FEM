/* $Id: FEManagerT_THK.cpp,v 1.13 2004-06-26 18:54:31 paklein Exp $ */
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
#include "iGridManagerT.h"
#include "iAutoArrayT.h"
#include "iNodeT.h"

using namespace Tahoe;

const double tol = 1.0e-3;   // for neighbor searching tolerance 1.0e-8 originally
const double root32 = sqrt(3.0)/2.0;    // for neighbor searching tolerance

/* constructor */
FEManagerT_THK::FEManagerT_THK(ifstreamT& input, ofstreamT& output, CommunicatorT& comm,
	const ArrayT<StringT>& argv, ifstreamT& bridging_input):
	FEManagerT_bridging(input, output, comm, argv, bridging_input)
{

}

/* initialize members */
void FEManagerT_THK::Initialize(InitCodeT init)
{
	/* inherited */
	FEManagerT_bridging::Initialize(init);
	
	ModelManagerT* model = FEManagerT::ModelManager();
	int nsd = model->NumDimensions();

	if (nsd == 2)
		Initialize2D();
	else if (nsd == 3)
		Initialize3D();
	else
	{
		const char caller[] = "FEManagerT_THK::Initialize";
		ExceptionT::GeneralFail(caller, "1D BRIDGING SCALE NOT ENABLED");
	}
	
}

/* 2D Bridging Scale Initialization */
void FEManagerT_THK::Initialize2D(void)
{
	ModelManagerT* model = FEManagerT::ModelManager();

	// read other parameters and initialize data
	ifstreamT& in = Input();
	in >> fNcrit;
	int num_neighbors = 2 * fNcrit + 1;   // maximum number of neighbors per atom in 2D
		  
	/* obtain list of atoms on which BC's will be applied in FEManagerT_THK */
	ArrayT<StringT> id_list;

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

	double lparam;
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
	ComputeThetaTables2D(data_file);
}

/* Bridging Scale 3D Initialization */
void FEManagerT_THK::Initialize3D(void)
{
	/* Implement 3D version of initialize here */
	/* Find neighbors of top 2 planes of atoms */
	ModelManagerT* model = FEManagerT::ModelManager();
	double lparam;
	int nsd = model->NumDimensions();

	// read other parameters and initialize data
	ifstreamT& in = Input();
	in >> fNcrit;
	int num_neighborsx = 2 * fNcrit + 1;	// assume same number of neighbors in both x and y directions
	int num_neighborsy = 2 * fNcrit + 1;
	fNeighbors = num_neighborsx * num_neighborsy;
	
	/* obtain list of atoms on which BC's will be applied in FEManagerT_THK */
	ArrayT<StringT> id_list;
        
	/* read node set indexes */
	model->NodeSetList(in, id_list);
	int numsets = id_list.Length();		// number of MD THK boundary node sets
 
	/* collect node sets - currently top/bottom 2 planes of ghost atoms only */
	fTopatoms = model->NodeSet(id_list[0]);	// row 3 ghost atoms (top)
	fTopatoms2 = model->NodeSet(id_list[1]);	// row 2 ghost atoms (top)
	fBottomatoms = model->NodeSet(id_list[2]);	// row 3 ghost atoms (bottom)
	fBottomatoms2 = model->NodeSet(id_list[3]);	// row 2 ghost atoms (bottom)
	
	/* read in top/bottom 2 planes of boundary atoms */
	fTA0 = model->NodeSet(id_list[4]);	// row 0 top atoms
	fTA1 = model->NodeSet(id_list[5]);	// row 1 top atoms
	fBA0 = model->NodeSet(id_list[6]);	// row 0 bottom atoms
	fBA1 = model->NodeSet(id_list[7]);	// row 1 bottom atoms
	fBottomrow1.Dimension(fBA1.Length());
	fBottomrow1.SetValueToPosition();
	fBottomrow0.Dimension(fBA0.Length());
	fBottomrow0 = fBottomrow1;
	fBottomrow0+=fBottomrow1.Length();
	fToprow0.Dimension(fTA0.Length());
	fToprow0 = fBottomrow0;
	fToprow0+=fBottomrow1.Length();
	fToprow1.Dimension(fTA1.Length());
	fToprow1=fToprow0;
	fToprow1+=fToprow0.Length();
	fTop20.Dimension(fTopatoms2.Length(), fNeighbors);
	fTop21.Dimension(fTopatoms2.Length(), fNeighbors);
	fTop30.Dimension(fTopatoms.Length(), fNeighbors);
	fTop31.Dimension(fTopatoms.Length(), fNeighbors);
	fBottom20.Dimension(fBottomatoms2.Length(), fNeighbors);
	fBottom21.Dimension(fBottomatoms2.Length(), fNeighbors);
	fBottom30.Dimension(fBottomatoms.Length(), fNeighbors);
	fBottom31.Dimension(fBottomatoms.Length(), fNeighbors);
	fTop20 = -1;  // -1 -> no neighbor
	fTop21 = -1;
	fTop30 = -1;
	fTop31 = -1;
	fBottom20 = -1;
	fBottom21 = -1;
	fBottom30 = -1;
	fBottom31 = -1;
	
	in >> lparam;
	
	NodeManagerT* node = FEManagerT::NodeManager();
	const dArray2DT& initcoords = node->InitialCoordinates();  // init coords for all atoms
	
	/* read in number of non-image atoms from file (real + ghost) */
	int niatoms;
	in >> niatoms;
	dArray2DT actualatoms(niatoms,3);
	iArrayT asdf(niatoms);
	asdf.SetValueToPosition();
	actualatoms.RowCollect(asdf, initcoords);
	//iArrayT asdf(initcoords.MajorDim());
	//asdf.SetValueToPosition();
	//dArray2DT actualatoms = initcoords;

	/* configure search grid - CURRENTLY SEARCHING ONLY NON-IMAGE ATOMS (REAL+GHOST) */
	iGridManagerT grid(10, 100, actualatoms, &asdf);
	grid.Reset();
	dArrayT acoord1, acoord2, acoord3, acoord4, ncoord1, ncoord2, ncoord3, x_atom;
	double xdist1, ydist1, mag, mag1, mag2, mag3, xmag, ymag;
	int blah2;
	
	/* bottom plane n=3 ghost atoms */
	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		initcoords.RowAlias(fBottomatoms[i], acoord1);
		/* candidate points for 3->0 planes */
		const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord1.Pointer(), 1.01*sqrt(10.0)*lparam);
	
		for (int j = 0; j < hitsa.Length(); j++)
		{			
			/* bottom distance 3->0 */
			ncoord1.Alias(nsd, hitsa[j].Coords());
			xmag = ncoord1[0]-acoord1[0];
			ymag = ncoord1[1]-acoord1[1];
			mag1 = sqrt(xmag*xmag+ymag*ymag);
						
			/* Get plane 0 neighbors for plane 3 bottom ghost atoms */
			if (abs(acoord1[2]-ncoord1[2]+3*lparam) < tol && mag1 < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-ncoord1[0]+acoord1[0]+lparam*k);
						ydist1 = abs(-ncoord1[1]+acoord1[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fBottom30(i,blah2) = hitsa[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}

	/* top plane n=3 ghost atoms */
	for (int i = 0; i < fTopatoms.Length(); i++)
	{
		initcoords.RowAlias(fTopatoms[i], acoord2);
		/* candidate points for 3->0 planes */
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(acoord2.Pointer(), 1.01*sqrt(10.0)*lparam);
	
		for (int j = 0; j < hits.Length(); j++)
		{
			/* top distance 3->0 */
			x_atom.Alias(nsd, hits[j].Coords());
			xmag = x_atom[0]-acoord2[0];
			ymag = x_atom[1]-acoord2[1];
			mag = sqrt(xmag*xmag+ymag*ymag);
			
			/* Get plane 0 neighbors for plane 3 top ghost atoms */
			if (abs(acoord2[2]-x_atom[2]-3*lparam) < tol && mag < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-x_atom[0]+acoord2[0]+lparam*k);
						ydist1 = abs(-x_atom[1]+acoord2[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fTop30(i,blah2) = hits[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}
	
	/* top plane n=3 ghost atoms */
	for (int i = 0; i < fTopatoms.Length(); i++)
	{
		initcoords.RowAlias(fTopatoms[i], acoord3);
		/* candidate points for 3->1 planes */
		const AutoArrayT<iNodeT>& hits2 = grid.HitsInRegion(acoord3.Pointer(), 1.01*sqrt(6.0)*lparam);
		for (int j = 0; j < hits2.Length(); j++)
		{
			/* top distance 3->1 */
			ncoord2.Alias(nsd, hits2[j].Coords());
			xmag = ncoord2[0]-acoord3[0];
			ymag = ncoord2[1]-acoord3[1];
			mag2 = sqrt(xmag*xmag+ymag*ymag);
		
			/* Get plane 1 neighbors for plane 3 top ghost atoms */
			if (abs(acoord3[2]-ncoord2[2]-2*lparam) < tol && mag2 < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-ncoord2[0]+acoord3[0]+lparam*k);
						ydist1 = abs(-ncoord2[1]+acoord3[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fTop31(i,blah2) = hits2[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}

	/* bottom plane n=3 ghost atoms */
	for (int i = 0; i < fBottomatoms.Length(); i++)
	{
		initcoords.RowAlias(fBottomatoms[i], acoord4);
		/* candidate points for 3->1 planes */
		const AutoArrayT<iNodeT>& hitsb = grid.HitsInRegion(acoord4.Pointer(), 1.01*sqrt(6.0)*lparam);
		for (int j = 0; j < hitsb.Length(); j++)
		{
			/* top distance 3->1 */
			ncoord3.Alias(nsd, hitsb[j].Coords());
			xmag = ncoord3[0]-acoord4[0];
			ymag = ncoord3[1]-acoord4[1];
			mag3 = sqrt(xmag*xmag+ymag*ymag);
					
			/* Get plane 1 neighbors for plane 3 bottom ghost atoms */
			if (abs(acoord4[2]-ncoord3[2]+2*lparam) < tol && mag3 < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-ncoord3[0]+acoord4[0]+lparam*k);
						ydist1 = abs(-ncoord3[1]+acoord4[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fBottom31(i,blah2) = hitsb[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}

	/* bottom plane n=2 ghost atoms */
	for (int i = 0; i < fBottomatoms2.Length(); i++)
	{
		initcoords.RowAlias(fBottomatoms2[i], acoord1);
		
		/* candidate points for 2->0 planes */
		const AutoArrayT<iNodeT>& hitsa = grid.HitsInRegion(acoord1.Pointer(), 1.01*sqrt(6.0)*lparam);
				
		for (int j = 0; j < hitsa.Length(); j++)
		{
			/* bottom distance 2->0 */
			ncoord1.Alias(nsd, hitsa[j].Coords());
			xmag = ncoord1[0]-acoord1[0];
			ymag = ncoord1[1]-acoord1[1];
			mag1 = sqrt(xmag*xmag+ymag*ymag);
						
			/* Get plane 0 neighbors for plane 2 bottom ghost atoms */
			if (abs(acoord1[2]-ncoord1[2]+2*lparam) < tol && mag1 < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-ncoord1[0]+acoord1[0]+lparam*k);
						ydist1 = abs(-ncoord1[1]+acoord1[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fBottom20(i,blah2) = hitsa[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}

	/* top plane n=2 ghost atoms */
	for (int i = 0; i < fTopatoms2.Length(); i++)
	{
		initcoords.RowAlias(fTopatoms2[i], acoord2);
		
		/* candidate points for 2->0 planes */
		const AutoArrayT<iNodeT>& hits = grid.HitsInRegion(acoord2.Pointer(), 1.01*sqrt(6.0)*lparam);
				
		for (int j = 0; j < hits.Length(); j++)
		{
			/* top distance 2->0 */
			x_atom.Alias(nsd, hits[j].Coords());
			xmag = x_atom[0]-acoord2[0];
			ymag = x_atom[1]-acoord2[1];
			mag = sqrt(xmag*xmag+ymag*ymag);
			
			/* Get plane 0 neighbors for plane 2 top ghost atoms */
			if (abs(acoord2[2]-x_atom[2]-2*lparam) < tol && mag < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-x_atom[0]+acoord2[0]+lparam*k);
						ydist1 = abs(-x_atom[1]+acoord2[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fTop20(i,blah2) = hits[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}
	
	/* top plane n=2 ghost atoms */
	for (int i = 0; i < fTopatoms2.Length(); i++)
	{
		initcoords.RowAlias(fTopatoms2[i], acoord3);
		/* candidate points for 2->1 planes */
		const AutoArrayT<iNodeT>& hits2 = grid.HitsInRegion(acoord3.Pointer(), 1.01*sqrt(2.0)*lparam);
		for (int j = 0; j < hits2.Length(); j++)
		{
			/* top distance 2->1 */
			ncoord2.Alias(nsd, hits2[j].Coords());
			xmag = ncoord2[0]-acoord3[0];
			ymag = ncoord2[1]-acoord3[1];
			mag2 = sqrt(xmag*xmag+ymag*ymag);

			/* Get plane 1 neighbors for plane 2 top ghost atoms */
			if (abs(acoord3[2]-ncoord2[2]-1*lparam) < tol && mag2 < sqrt(2.0)*1.01*lparam)	
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-ncoord2[0]+acoord3[0]+lparam*k);
						ydist1 = abs(-ncoord2[1]+acoord3[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fTop21(i,blah2) = hits2[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}
	
	/* bottom plane n=2 ghost atoms */
	for (int i = 0; i < fBottomatoms2.Length(); i++)
	{
		initcoords.RowAlias(fBottomatoms2[i], acoord4);
		/* candidate points for 2->1 planes */
		const AutoArrayT<iNodeT>& hits2 = grid.HitsInRegion(acoord4.Pointer(), 1.01*sqrt(2.0)*lparam);
		for (int j = 0; j < hits2.Length(); j++)
		{
			/* bottom distance 2->1 */
			ncoord3.Alias(nsd, hits2[j].Coords());
			xmag = ncoord3[0]-acoord4[0];
			ymag = ncoord3[1]-acoord4[1];
			mag2 = sqrt(xmag*xmag+ymag*ymag);

			/* Get plane 1 neighbors for plane 2 bottom ghost atoms */
			if (abs(acoord4[2]-ncoord3[2]+1*lparam) < tol && mag3 < sqrt(2.0)*1.01*lparam)
			{
				/* now sort into x by y array for THK BC application */
				blah2 = 0;
				for (int k = -fNcrit; k <= fNcrit; k++)
				{
					for (int l = -fNcrit; l<= fNcrit; l++)
					{
						xdist1 = abs(-ncoord3[0]+acoord4[0]+lparam*k);
						ydist1 = abs(-ncoord3[1]+acoord4[1]+lparam*l);
						if (xdist1 < tol && ydist1 < tol)
							fBottom21(i,blah2) = hits2[j].Tag();
						blah2++;
					}	
				}
			}
		}
	}
	
	/* Compute Theta tables - may need to pass multiple data files as arguments */
	StringT path, path2, path3, path4;
	path.FilePath(in.filename());
	StringT data_file, data_file2, data_file3, data_file4;
	in >> data_file;
	data_file.ToNativePathName();
	data_file.Prepend(path);
	path2.FilePath(in.filename());
	in >> data_file2;
	data_file2.ToNativePathName();
	data_file2.Prepend(path2);	
	in >> data_file3;
	data_file3.ToNativePathName();
	data_file3.Prepend(path3);	
	in >> data_file4;
	data_file4.ToNativePathName();
	data_file4.Prepend(path4);	
	ComputeThetaTables3D(data_file, data_file2, data_file3, data_file4);
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

/* return iArrayT of boundary and ghost atom numbers - 3D version */
const iArrayT& FEManagerT_THK::InterpolationNodes3D(void)
{
	const iArrayT& ghost = GhostNodes();
	int size = ghost.Length()+fBA0.Length()+fBA1.Length()+fTA0.Length()+fTA1.Length();
	fInterpolationNodes.Dimension(size);
	fInterpolationNodes.CopyIn(0,ghost);	// copy in ghost atoms first
	fInterpolationNodes.CopyIn(ghost.Length(), fBA1);	// then copy in bottom atoms row 1
	fInterpolationNodes.CopyIn(ghost.Length()+fBA1.Length(), fBA0);	// bottom atoms row 0
	fInterpolationNodes.CopyIn(ghost.Length()+fBA1.Length()+fBA0.Length(), fTA0);	// top atoms row 0
	fInterpolationNodes.CopyIn(ghost.Length()+fBA1.Length()+fBA0.Length()+fTA0.Length(), fTA1);	// top atoms row 1
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

/* calculate ghost atom fine scale displacements using 3D disp/disp formulation */
const dArray2DT& FEManagerT_THK::THKDisp(const dArray2DT& badisp)
{
	StringT bridging_field = "displacement";
	
	/* badisp is in format bottom displacements first, then top */
	/* COULD DIMENSION ONCE IN INTERPOLATIONNODES3D() */
	fTHKdisp.Dimension(badisp.MajorDim(), 3);  
	fTHKdisp = 0.0;

	/* sort top and bottom displacements into separate arrays - may not be necessary */
	/* ADD topdisp0, etc as class variables so only dimension once? */
	/* COULD DIMENSION ONCE IN INTERPOLATIONNODES3D() */
	dArray2DT topdisp0, topdisp1, bottomdisp0, bottomdisp1;
	dArrayT atomdisp, femdisp(3), diff(3);
	bottomdisp0.Dimension(fBA0.Length(), 3);
	bottomdisp1.Dimension(fBA1.Length(), 3);
	topdisp0.Dimension(fTA0.Length(), 3);
	topdisp1.Dimension(fTA1.Length(), 3);
	bottomdisp0.RowCollect(fBottomrow0, badisp);
	bottomdisp1.RowCollect(fBottomrow1, badisp);
	topdisp0.RowCollect(fToprow0, badisp);
	topdisp1.RowCollect(fToprow1, badisp);
	
	/* access the actual MD displacements */
	NodeManagerT* node = FEManagerT::NodeManager();
	FieldT* atomfield = node->Field(bridging_field);
	dArray2DT mddisp = (*atomfield)[0];

	const int stepnum = FEManagerT::StepNumber();  // to write into correct part of fHistoryTable 
	const double timestep = FEManagerT::TimeStep();  // delta t_md

	dArray2DT shift;
	shift.Dimension(fNumstep_crit-1, 3);	// copy all rows of history except last row 

	/* Calculate q - ubar for top/bottom plane 0 of atoms */
	for (int i = 0; i < fTA0.Length(); i++)
	{
		/* non-shift case */
		if (stepnum < fNumstep_crit)   
		{
			/* Row 0 top plane of atoms */
			mddisp.RowAlias(fTA0[i], atomdisp);
			topdisp0.RowAlias(i, femdisp);	
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTablet[i].SetRow(stepnum, diff);
		
			/* Row 0 bottom plane of atoms */
			mddisp.RowAlias(fBA0[i], atomdisp);
			bottomdisp0.RowAlias(i, femdisp);	
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTableb[i].SetRow(stepnum, diff);
		}
		else	// t > t_crit
		{
			/* Row 0 top plane of atoms */
			mddisp.RowAlias(fTA0[i], atomdisp);
			topdisp0.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTablet[i]);
			fHistoryTablet[i].BlockRowCopyAt(shift, 0);
			fHistoryTablet[i].SetRow(fNumstep_crit-1, diff);
	
			/* Row 0 bottom plane of atoms */
			mddisp.RowAlias(fBA0[i], atomdisp);
			bottomdisp0.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTableb[i]);
			fHistoryTableb[i].BlockRowCopyAt(shift, 0);
			fHistoryTableb[i].SetRow(fNumstep_crit-1, diff);
		}
	}

	/* Calculate q - ubar for top/bottom plane 1 of atoms */
	for (int i = 0; i < fTA1.Length(); i++)
	{
		/* non-shift case */
		if (stepnum < fNumstep_crit)   
		{
			/* Row 1 top plane of atoms */
			mddisp.RowAlias(fTA1[i], atomdisp);
			topdisp1.RowAlias(i, femdisp);	
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTablet1[i].SetRow(stepnum, diff);
		
			/* Row 1 bottom plane of atoms */
			mddisp.RowAlias(fBA1[i], atomdisp);
			bottomdisp1.RowAlias(i, femdisp);	
			diff.DiffOf(atomdisp, femdisp);
			fHistoryTableb1[i].SetRow(stepnum, diff);
		}
		else	// t > t_crit
		{
			/* Row 1 top plane of atoms */
			mddisp.RowAlias(fTA1[i], atomdisp);
			topdisp1.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTablet1[i]);
			fHistoryTablet1[i].BlockRowCopyAt(shift, 0);
			fHistoryTablet1[i].SetRow(fNumstep_crit-1, diff);
	
			/* Row 1 bottom plane of atoms */
			mddisp.RowAlias(fBA1[i], atomdisp);
			bottomdisp1.RowAlias(i, femdisp);
			diff.DiffOf(atomdisp, femdisp);
			shift.RowCollect(fShift, fHistoryTableb1[i]);
			fHistoryTableb1[i].BlockRowCopyAt(shift, 0);
			fHistoryTableb1[i].SetRow(fNumstep_crit-1, diff);
		}
	}

	dMatrixT theta(3);
	dArrayT force1(3), force2(3), force0a(3), force0b(3), force0c(3), force0d(3);
	force0c = 0.0;
	force0d = 0.0;
	force0a = 0.0;
	force0b = 0.0;
	
	/* set inverse maps */
	InverseMapT topnodes0, bottomnodes0, topnodes1, bottomnodes1;
	topnodes1.SetMap(fTA1);	// Set global to local map for top atoms plane 1
	bottomnodes1.SetMap(fBA1);	// Set global to local map for bottom atoms plane 1
	topnodes0.SetMap(fTA0);	// Set global to local map for top atoms plane 0
	bottomnodes0.SetMap(fBA0);	// Set global to local map for bottom atoms plane 0
	int count, dex2, dex, dex3, dex4;

	/* loop over ghost atom neighbors, calculate THK disp, i.e. Theta(t-tau)*(q-ubar)(tau)*(delta tau) */
	/* calculate THK disp for plane 2 top/bottom ghost atoms */
	for (int i = 0; i < fTop20.MajorDim(); i++)
	{
		count = 0;
		for (int j = 0; j < 2*fNcrit+1; j++)
		{
			for (int k = 0; k < 2*fNcrit+1; k++)
			{
				if (fTop20(i,count) == -1)	// assume fTop20=fBottom20
					int blah = 0;
				else
				{
					/* top plane 0 */
					dex2 = topnodes0.Map(fTop20(i,count));  
					const dArray2DT& disp_t0 = fHistoryTablet[dex2];
					const dArray2DT& theta_t0 = fTheta11t[fNeighbors-1-count]; 	// start count at 8, go backwards
					
					/* bottom plane 0 */
					dex = bottomnodes0.Map(fBottom20(i,count));
					const dArray2DT& disp_b0 = fHistoryTableb[dex];
					const dArray2DT& theta_b0 = fTheta11b[fNeighbors-1-count]; 	// start count at 8, go backwards 
					
					/* calculate fine scale THK disp using theta here */
					if (stepnum < fNumstep_crit)  // < fNumstep_crit 
					{
						for (int l = 0; l < stepnum; l++)	// ORIGINALLY l < stepnum
						{
							/* top plane 0 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0a.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{11t}*u_{l',m',0}
												
							/* bottom plane 0 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{11b}*u_{l',m',0}
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{
							/* top plane 0 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0a.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{11t}*u_{l',m',0}
												
							/* bottom plane 0 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{11b}*u_{l',m',0}
						}
					}
				}
				
				if (fTop21(i,count) == -1)	// assume fTop21=fBottom21
					int blah = 0;
				else
				{
					/* top plane 0 */
					dex2 = topnodes1.Map(fTop21(i,count));  
					const dArray2DT& disp_t0 = fHistoryTablet1[dex2];
					const dArray2DT& theta_t0 = fTheta12t[fNeighbors-1-count]; 	// start count at 8, go backwards
		
					/* bottom plane 0 */
					dex = bottomnodes1.Map(fBottom21(i,count));  
					const dArray2DT& disp_b0 = fHistoryTableb1[dex];
					const dArray2DT& theta_b0 = fTheta12b[fNeighbors-1-count]; 	// start count at 8, go backwards 
				
					/* calculate fine scale THK disp using theta here */
					if (stepnum < fNumstep_crit)  // < fNumstep_crit 
					{
						for (int l = 0; l < stepnum; l++)	// ORIGINALLY l < stepnum
						{
							/* top plane 0 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0c.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{12t}*u_{l',m',1}
												
							/* bottom plane 0 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0d.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{12b}*u_{l',m',1}
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{
							/* top plane 0 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0c.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{12t}*u_{l',m',1}
											
							/* bottom plane 0 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0d.AddScaled(timestep, force2);	// u_{l,m,2}=Theta^{12b}*u_{l',m',1}
						}
					}
				}
				count++;
			}
		}

		/* add THK GA disp to fTHKdisp - top and bottom atoms 0 */
		fTHKdisp.AddToRowScaled(fBottomrow0[i], 1.0, force0b);
		fTHKdisp.AddToRowScaled(fToprow0[i], 1.0, force0a);
		
		/* add THK GA disp to fTHKdisp - top and bottom atoms 1 */
		fTHKdisp.AddToRowScaled(fBottomrow0[i], 1.0, force0d);
		fTHKdisp.AddToRowScaled(fToprow0[i], 1.0, force0c);
	}

	/* NEED TO ZERO OUT FORCE0a/b/c/d here? */
	force0a = 0.0;
	force0b = 0.0;
	force0c = 0.0;
	force0d = 0.0;

	/* loop over ghost atom neighbors, calculate THK disp, i.e. Theta(t-tau)*(q-ubar)(tau)*(delta tau) */
	/* calculate THK disp for plane 3 top/bottom ghost atoms */
	for (int i = 0; i < fTop31.MajorDim(); i++)
	{
		count = 0;
		for (int j = 0; j < 2*fNcrit+1; j++)
		{
			for (int k = 0; k < 2*fNcrit+1; k++)
			{
				if (fTop30(i,count) == -1)	// currently assuming fTop30 = fBottom30 in terms of neighbors
					int blah = 0;
				else
				{
					/* top plane 1 */
					dex2 = topnodes0.Map(fTop30(i,count));  
					const dArray2DT& disp_t0 = fHistoryTablet[dex2];
					const dArray2DT& theta_t0 = fTheta21t[fNeighbors-1-count];
					
					/* bottom plane 1 */
					dex = bottomnodes0.Map(fBottom30(i,count));  
					const dArray2DT& disp_b0 = fHistoryTableb[dex];
					const dArray2DT& theta_b0 = fTheta21b[fNeighbors-1-count];
					
					/* calculate fine scale THK disp using theta here */
					if (stepnum < fNumstep_crit)  // < fNumstep_crit 
					{
						for (int l = 0; l < stepnum; l++)	// ORIGINALLY l < stepnum
						{
							/* top plane 1 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0a.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{21t}*u_{l',m',0}
												
							/* bottom plane 1 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{21b}*u_{l',m',0}
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{
							/* top plane 1 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0a.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{21t}*u_{l',m',0}
												
							/* bottom plane 1 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0b.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{21b}*u_{l',m',0}
						}
					}
				}
				if (fTop31(i,count) == -1)	// currently assuming fTop31 = fBottom31 in terms of neighbors
					int blah = 0;
				else
				{
					/* top plane 1 */
					dex2 = topnodes1.Map(fTop31(i,count));  
					const dArray2DT& disp_t0 = fHistoryTablet1[dex2];
					const dArray2DT& theta_t0 = fTheta22t[fNeighbors-1-count];
					
					/* bottom plane 1 */
					dex = bottomnodes1.Map(fBottom31(i,count));  
					const dArray2DT& disp_b0 = fHistoryTableb1[dex];
					const dArray2DT& theta_b0 = fTheta22b[fNeighbors-1-count];
					
					/* calculate fine scale THK disp using theta here */
					if (stepnum < fNumstep_crit)  // < fNumstep_crit 
					{
						for (int l = 0; l < stepnum; l++)	// ORIGINALLY l < stepnum
						{
							/* top plane 1 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0c.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{22t}*u_{l',m',1}
												
							/* bottom plane 1 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(stepnum-l, force1);
							theta.Multx(force1, force2);
							force0d.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{22b}*u_{l',m',1}
						}
					}
					else	// normalized time greater than critical value
					{	
						for (int l = 0; l < fNumstep_crit; l++)
						{
							/* top plane 1 */
							theta.Alias(3,3,theta_t0(l));	// MAY NEED TO USE l+1 here 
							disp_t0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0c.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{22t}*u_{l',m',1}
												
							/* bottom plane 1 */
							theta.Alias(3,3,theta_b0(l));	// MAY NEED TO USE l+1 here 
							disp_b0.RowAlias(fNumstep_crit-l-1, force1);
							theta.Multx(force1, force2);
							force0d.AddScaled(timestep, force2);	// u_{l,m,3}=Theta^{22b}*u_{l',m',1}
						}
					}
				}
				count++;
			}
		}

		/* add THK GA disp to fTHKdisp - top and bottom atoms 0 */
		fTHKdisp.AddToRowScaled(fBottomrow1[i], 1.0, force0b);
		fTHKdisp.AddToRowScaled(fToprow1[i], 1.0, force0a);
		
		/* add THK GA disp to fTHKdisp - top and bottom atoms 1 */
		fTHKdisp.AddToRowScaled(fBottomrow1[i], 1.0, force0d);
		fTHKdisp.AddToRowScaled(fToprow1[i], 1.0, force0c);
	}

	//cout << "THKdisp = " << fTHKdisp << endl;
	return fTHKdisp;

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

/* function to modify initial positions via Gaussian displacement in y-direction */
const dArray2DT& FEManagerT_THK::ThreeDWave(const iArrayT& nonghostatoms) 
{
	double sig = 15.0;
	double H = sig/4.0;
	double A = 1.5e-2;
	double b = 0.1;
	double Lc = 4.0 * sig;
	double uc = exp(-pow(Lc/sig,2));
	double pi = acos(-1.0);
	double xi, yi, zi, dx, dy, dz, dist;
	
	/* loop over all atoms, apply initial displacement */
	NodeManagerT* node = FEManagerT::NodeManager();
	const dArray2DT& initcoords = node->InitialCoordinates();  // init coords for all atoms
	dArray2DT initdisp, coord1;
	coord1.Dimension(nonghostatoms.Length(), initcoords.MinorDim());
	fInitdisp.Dimension(nonghostatoms.Length(), initcoords.MinorDim());
	coord1.RowCollect(nonghostatoms, initcoords);
	
	for (int i = 0; i < nonghostatoms.Length(); i++)
	{
		xi = fabs(coord1(i,0));
		yi = fabs(coord1(i,1));
		zi = fabs(coord1(i,2));
		dist = sqrt(xi*xi+yi*yi+zi*zi);	
		if (dist <= Lc)
		{
			dx = A/(1.0-uc)*(1.0+b*cos(2*pi*xi/H))*(exp(-pow(xi/sig,2))-uc);
			dy = A/(1.0-uc)*(1.0+b*cos(2*pi*yi/H))*(exp(-pow(yi/sig,2))-uc);
			dz = A/(1.0-uc)*(1.0+b*cos(2*pi*zi/H))*(exp(-pow(zi/sig,2))-uc);
			fInitdisp(i,0) = dx;
			fInitdisp(i,1) = dy;
			fInitdisp(i,2) = dz;
		}
	}
	return fInitdisp;
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
void FEManagerT_THK::ComputeThetaTables3D(const StringT& data_file, const StringT& data_file2, const StringT& data_file3, const StringT& data_file4)
{
	const char caller[] = "FEManagerT_THK::ComputeThetaTables3D";
	ifstreamT data(data_file);
	if (!data.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file.Pointer());
		
	ifstreamT data2(data_file2);
	if (!data2.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file2.Pointer());
	
	ifstreamT data3(data_file3);
	if (!data3.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file3.Pointer());
		
	ifstreamT data4(data_file4);
	if (!data4.is_open())
		ExceptionT::GeneralFail(caller, "file not found: %s", data_file4.Pointer());
				
	/* dimensions */
	double pi = acos(-1.0);
	int n_neighbor = (2 * fNcrit + 1) * (2 * fNcrit + 1);  // number of neighbors to calculate thetas for
	int n_sum, nsteps;       // number of fourier coefficients to use in calculation of theta
	data >> n_sum;
	data2 >> n_sum;
	data3 >> n_sum;
	data4 >> n_sum;
	
	/* dimension work space */
	ifstreamT& in = Input();
	double tstep, totaltime, looptime;  // timestep
	in >> fN_times;
	in >> tstep;
	totaltime = fN_times * tstep;   // total time of simulation
	
	/* determine correct loop time for theta and time history variables */
	if (totaltime <= 12.0)  // assume that t_crit = 12.0 in normalized time
		looptime = totaltime;
	else
		looptime = 12.0;   // need to normalize this time
	fNumstep_crit = int(12.0/tstep) + 1;	// currently assuming normalized t_crit = 12.0 
	
	/* determine correct number of timesteps to store for theta and history variables */
	if (fN_times+1 <= fNumstep_crit)	// WATCH +1 FACTOR!!!
		nsteps = fN_times+1;  // watch +1
	else
		nsteps = fNumstep_crit;

	iArrayT crit(fNumstep_crit);
	fShift.Dimension(fNumstep_crit-1);
	crit.SetValueToPosition();
	fShift.CopyPart(0, crit, 1, fNumstep_crit-1);	// Used to shift displacement history when time > t_crit - MAY NEED TO FIX!
	
	fTheta11t.Dimension(n_neighbor);
	fTheta11b.Dimension(n_neighbor);
	fTheta12t.Dimension(n_neighbor);
	fTheta12b.Dimension(n_neighbor);	// may need more ThetaTables (11, 12, 21, 22) for both top/bottom
	fTheta21t.Dimension(n_neighbor);
	fTheta21b.Dimension(n_neighbor);
	fTheta22t.Dimension(n_neighbor);
	fTheta22b.Dimension(n_neighbor);
	
	/* dimension boundary atom history */
	fHistoryTablet.Dimension(fTA0.Length());
	fHistoryTableb.Dimension(fBA0.Length());
	fHistoryTablet1.Dimension(fTA1.Length());
	fHistoryTableb1.Dimension(fBA1.Length());	// Do not need more boundary atom histories
	
	for (int i = 0; i < fTA0.Length(); i++)
	{
		fHistoryTablet[i].Dimension(nsteps, 3);  // used to be fN_times+1,2
		fHistoryTableb[i].Dimension(nsteps, 3);
	}
	
	for (int i = 0; i < fTA1.Length(); i++)
	{
		fHistoryTablet1[i].Dimension(nsteps, 3);  // used to be fN_times+1,2
		fHistoryTableb1[i].Dimension(nsteps, 3);
	}
	
	/* dimension theta table data structure */
	for (int i = 0; i < fTheta11t.Length(); i++)
	{
		fTheta11t[i].Dimension(nsteps, 3*3);  // used to be fN_times+1,3*3
		fTheta11b[i].Dimension(nsteps, 3*3);
		fTheta12t[i].Dimension(nsteps, 3*3);
		fTheta12b[i].Dimension(nsteps, 3*3);
		fTheta21t[i].Dimension(nsteps, 3*3);
		fTheta21b[i].Dimension(nsteps, 3*3);
		fTheta22t[i].Dimension(nsteps, 3*3);
		fTheta22b[i].Dimension(nsteps, 3*3);
	}
	
	/* read in data tables for multiple thetas */
	dArrayT row, row2, row3, row4;
	ArrayT<dArray2DT> data_table(n_neighbor), data_table2(n_neighbor), data_table3(n_neighbor), data_table4(n_neighbor);
	for (int i = 0; i < data_table.Length(); i++)
	{
		/* dimension */
		dArray2DT& n_table = data_table[i];
		dArray2DT& n_table2 = data_table2[i];
		dArray2DT& n_table3 = data_table3[i];
		dArray2DT& n_table4 = data_table4[i];
		n_table.Dimension(n_sum, 3*3);
		n_table2.Dimension(n_sum, 3*3);
		n_table3.Dimension(n_sum, 3*3);
		n_table4.Dimension(n_sum, 3*3);
				
		/* read */
		for (int j = 0; j < n_table.MajorDim(); j++)
		{
			int junk, junk2, junk3, junk4;
			data >> junk >> junk >> junk;
			data2 >> junk2 >> junk2 >> junk2;
			data3 >> junk3 >> junk3 >> junk3;
			data4 >> junk4 >> junk4 >> junk4;
			
			/* each row is b^T */
			n_table.RowAlias(j, row);
			n_table2.RowAlias(j, row2);
			n_table3.RowAlias(j, row3);
			n_table4.RowAlias(j, row4);
			data >> row;
			data2 >> row2;
			data3 >> row3;
			data4 >> row4;
		}
	}
		
	dMatrixT theta1(3), temptheta(3), temptheta2(3), temptheta3(3), temptheta4(3), transform(3);
	transform.Identity(1.0);
	transform(2,2) = -1.0; // Theta(bottom)=T*Theta(top)*T, T=transform

	/* compute theta's - top/bottom planes related by transform matrix above */
	for (int i = 0; i < n_neighbor; i++)	// fThetaTable.Length()
	{
		int count = 0;
		/* get each set of Fourier coefficients */
		const dArray2DT& theta_i = data_table[i];
		const dArray2DT& theta_j = data_table2[i];
		const dArray2DT& theta_k = data_table3[i];
		const dArray2DT& theta_l = data_table4[i];
		
		for (double j = 0.0; j <= looptime; j+=tstep)	// theta only goes to normalized time = 3
		{
			int n_theta = data_table[i].MajorDim();
			temptheta = 0.0;
			temptheta2 = 0.0;
			temptheta3 = 0.0;
			temptheta4 = 0.0;
			
			for (int k = 0; k < n_theta; k++)  // can truncate this summation
			{
				/* Extract each row of coefficients */
				temptheta.AddScaled(sin((k+1)*pi*j/12.2), theta_i(k)); // tmax = 12.2 now
				temptheta2.AddScaled(sin((k+1)*pi*j/12.2), theta_j(k)); // tmax = 12.2 now
				temptheta3.AddScaled(sin((k+1)*pi*j/12.2), theta_k(k)); // tmax = 12.2 now
				temptheta4.AddScaled(sin((k+1)*pi*j/12.2), theta_l(k)); // tmax = 12.2 now
			}

			/* add temptheta into fThetaTable */
			fTheta11t[i].SetRow(count, temptheta);
			fTheta12t[i].SetRow(count, temptheta2);
			fTheta21t[i].SetRow(count, temptheta3);
			fTheta22t[i].SetRow(count, temptheta4);
		
			/* compute bottom thetatables via transform matrix */
			theta1.MultABC(transform, temptheta, transform);
			fTheta11b[i].SetRow(count, theta1);
			theta1.MultABC(transform, temptheta2, transform);
			fTheta12b[i].SetRow(count, theta1);
			theta1.MultABC(transform, temptheta3, transform);
			fTheta21b[i].SetRow(count, theta1);
			theta1.MultABC(transform, temptheta4, transform);
			fTheta22b[i].SetRow(count, theta1);
			count++;
		}
	}

	//for (int i = 0; i < n_neighbor; i++)
	//{
	//	cout << "Thetatable11 TOP = " << fTheta11t[i] << endl;
	//	cout << "Thetatable11 BOTTOM = " << fTheta11b[i] << endl;
	//	cout << "Thetatable12 TOP = " << fTheta12t[i] << endl;
	//	cout << "Thetatable12 BOTTOM = " << fTheta12b[i] << endl;
	//	cout << "Thetatable21 TOP = " << fTheta21t[i] << endl;
	//	cout << "Thetatable21 BOTTOM = " << fTheta21b[i] << endl;
	//	cout << "Thetatable22 TOP = " << fTheta22t[i] << endl;
	//	cout << "Thetatable22 BOTTOM = " << fTheta22b[i] << endl;
	//}
}

#endif /* BRIDGING_ELEMENT */
