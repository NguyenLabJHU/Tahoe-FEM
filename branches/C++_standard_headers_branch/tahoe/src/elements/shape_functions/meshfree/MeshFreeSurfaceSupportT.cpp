/* $Id: MeshFreeSurfaceSupportT.cpp,v 1.5 2003-11-21 22:47:14 paklein Exp $ */
/* created: paklein (02/22/2000)                                          */
/* supporting functions for cohesive elements in a meshfree domain        */

#include "MeshFreeSurfaceSupportT.h"

#include "ExceptionT.h"
#include "toolboxConstants.h"

#include "SurfaceShapeT.h"
#include "MeshFreeSupportT.h"
#include "LocalArrayT.h"

//DEBUG
//#include <fstream.h>
//DEBUG

/* parameters */

using namespace Tahoe;

const int kHeadRoom = 20; // percent

/* constructor */
MeshFreeSurfaceSupportT::MeshFreeSurfaceSupportT(MeshFreeSupportT& mf_support,
	SurfaceShapeT& ref_surface_shape, LocalArrayT& ref_loc_coords,
	const dArray2DT& facet_coords, int num_facet_nodes, bool storeshape):
	fStoreShape(storeshape),
	fMFSupport(mf_support),
	fRefSurfaceShape(ref_surface_shape),
	fRefLocCoords(ref_loc_coords),
	fNumFacetNodes(num_facet_nodes),
	fFacetCoords(facet_coords),
	fNeighbors_1(kHeadRoom), fPhi_1(kHeadRoom), fDPhi_1(kHeadRoom),
	fNeighbors_2(kHeadRoom), fPhi_2(kHeadRoom), fDPhi_2(kHeadRoom)
{
//TEMP - current LoadData requires stored shapefunctions
	if (!fStoreShape)
	{
		cout << "\n MeshFreeSurfaceSupportT::MeshFreeSurfaceSupportT: setting\n"
		     <<   "     storage option to TRUE" << endl;
		fStoreShape = true;
	}
}

/* (re-)compute shape functions - NULL resets all facets */
void MeshFreeSurfaceSupportT::ResetFacets(const ArrayT<int>* reset_facets)
{
	/* with storage only */
	if (!fStoreShape) return;	

	/* check */
	if (fRefLocCoords.NumberOfNodes() !=
	    fFacetCoords.MinorDim()/(fRefSurfaceShape.NumSD() + 1))
	{
		cout << "\n MeshFreeSurfaceSupportT::ResetFacets: size mismatch with local\n"
		     <<   "     facet coords." << endl;
		throw ExceptionT::kSizeMismatch;
	}

	/* dimensions */
	int nip = fRefSurfaceShape.NumIP();
	int nsd = fRefSurfaceShape.NumSD() + 1;
	int nft = fFacetCoords.MajorDim();

	/* storage for ip data */	
	ArrayT< AutoArrayT<int> > neighbors_1(nip);
	ArrayT< AutoArrayT<double> > phi_1(nip);
	ArrayT< AutoArrayT<double> > Dphi_1(nip);

	ArrayT< AutoArrayT<int> > neighbors_2(nip);
	ArrayT< AutoArrayT<double> > phi_2(nip);
	ArrayT< AutoArrayT<double> > Dphi_2(nip);
	
	dArray2DT facet_coords;
	dMatrixT Q(nsd);
	dArrayT normal(nsd);

//DEBUG
#if 0
const dArray2DT& coordinates = fMFSupport.NodalCoordinates();
ofstream out1("side1.out");
ofstream out2("side2.out");
out1 << " MeshFreeSurfaceSupportT::SetShapeFunctions\n " << '\n';
out2 << " MeshFreeSurfaceSupportT::SetShapeFunctions\n " << '\n';
#endif
//DEBUG

	/* selectively or all */
	int lim = (reset_facets != NULL) ? reset_facets->Length() : nft;
	for (int j = 0; j < lim; j++)
	{
		int facet = (reset_facets != NULL) ? (*reset_facets)[j] : j;

//DEBUG
#if 0
out1 << "\n facet = " << facet+1 << '\n';
out2 << "\n facet = " << facet+1 << '\n';
#endif
//DEBUG

		/* set facet coordinates */
		facet_coords.Alias(fRefLocCoords.NumberOfNodes(), nsd, fFacetCoords(facet));
		fRefLocCoords.FromTranspose(facet_coords);
		double h = FacetSize(facet_coords);		

		/* loop over integration points */
		fRefSurfaceShape.TopIP();
		while (fRefSurfaceShape.NextIP())
		{
			int ip = fRefSurfaceShape.CurrIP();

//DEBUG
#if 0
out1 << "\n ip = " << ip+1 << '\n';
out2 << "\n ip = " << ip+1 << '\n';
#endif
//DEBUG
		
			/* (reference) coordinate transformation */
			fRefSurfaceShape.Jacobian(Q);

			/* ip coordinates */
			const dArrayT& ip_coord = fRefSurfaceShape.IPCoords();

			/* side 1: (-) */
			Q.CopyColumn(nsd-1, normal);
			normal *= -h/25.0;
			if (!fMFSupport.SetFieldAt(ip_coord, &normal))
			{
				cout << "\n MeshFreeSurfaceSupportT::SetShapeFunctions: could not form MLS on -:\n";
				cout <<   "     facet: " << facet+1 << '\n';
				cout <<   "        ip: " << ip+1 << '\n';
				cout <<   "         x: " << ip_coord.no_wrap() << '\n';
				cout <<   "         n: " << normal.no_wrap() << '\n';
				throw ExceptionT::kGeneralFail;
			}

			/* collect data */
			neighbors_1[ip] = fMFSupport.NeighborsAt();
			phi_1[ip]  = fMFSupport.FieldAt();
			Dphi_1[ip] = fMFSupport.DFieldAt();

//DEBUG
#if 0
for (int i = 0; i < neighbors.Length(); i++)
	out1 << setw(kIntWidth)    << neighbors[i]+1
	     << setw(kDoubleWidth) << (phi_1[ip])[i]
	     << setw(kDoubleWidth) << coordinates(neighbors[i],2) << '\n';
#endif
//DEBUG

			/* side 2: (+) */
			Q.CopyColumn(nsd-1, normal);
			normal *= h/25.0;
			if (!fMFSupport.SetFieldAt(ip_coord, &normal))
			{
				cout << "\n MeshFreeSurfaceSupportT::SetShapeFunctions: could not form MLS on +:\n";
				cout <<   "     facet: " << facet+1 << '\n';
				cout <<   "        ip: " << ip+1 << '\n';
				cout <<   "         x: " << ip_coord.no_wrap() << '\n';
				cout <<   "         n: " << normal.no_wrap() << '\n';
				throw ExceptionT::kGeneralFail;
			}

			/* collect data */
			neighbors_2[ip] = fMFSupport.NeighborsAt();
			phi_2[ip]  = fMFSupport.FieldAt();
			Dphi_2[ip] = fMFSupport.DFieldAt();

//DEBUG
#if 0
for (int i = 0; i < neighbors.Length(); i++)
	out2 << setw(kIntWidth)    << neighbors[i]+1
	     << setw(kDoubleWidth) << (phi_2[ip])[i]
	     << setw(kDoubleWidth) << coordinates(neighbors[i],2) << '\n';

#endif
//DEBUG
		}

		/* write to databases */
		StoreData(facet, 0, neighbors_1, phi_1, Dphi_1);
		StoreData(facet, 1, neighbors_2, phi_2, Dphi_2);
	}
}

void MeshFreeSurfaceSupportT::LoadData(int facet, int side, iArrayT& neighbors,
	dArray2DT& phi, ArrayT<dArray2DT>& Dphi)
{
#if __option(extended_errorcheck)
	if (Dphi.Length() != fRefSurfaceShape.NumIP()) throw ExceptionT::kSizeMismatch;
	if (side != 0 && side != 1) throw ExceptionT::kOutOfRange;
#endif

	/* dimensions */
	int nip = fRefSurfaceShape.NumIP();
	int nsd = fRefSurfaceShape.NumSD() + 1;

	if (fStoreShape)
	{
		/* neighbors */
		if (side == 0)
			fNeighbors_1.RowAlias(facet, neighbors);
		else
			fNeighbors_2.RowAlias(facet, neighbors);
		int nnd = neighbors.Length();
	
		/* load functions */
		phi.Set(nip, nnd, (side == 0) ? fPhi_1(facet) : fPhi_2(facet));
			
		/* load derivatives */
		double* Dphi_ptr = (side == 0) ? fDPhi_1(facet) : fDPhi_2(facet);
		for (int i = 0; i < nip; i++)
		{
			Dphi[i].Set(nsd, nnd, Dphi_ptr);
			Dphi_ptr += nsd*nnd;
		}
	}
	else
	{
		cout << "\n MeshFreeSurfaceSupportT::LoadData: must use stored functions" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

void MeshFreeSurfaceSupportT::LoadData(int facet,
	iArrayT& neighbors_1, dArray2DT& phi_1, ArrayT<dArray2DT>& Dphi_1,
	iArrayT& neighbors_2, dArray2DT& phi_2, ArrayT<dArray2DT>& Dphi_2)
{
#if __option(extended_errorcheck)
	if (Dphi_1.Length() != fRefSurfaceShape.NumIP() ||
		Dphi_2.Length() != fRefSurfaceShape.NumIP()) throw ExceptionT::kSizeMismatch;
#endif

	/* dimensions */
	int nip = fRefSurfaceShape.NumIP();
	int nsd = fRefSurfaceShape.NumSD() + 1;

	if (fStoreShape)
	{
		/* neighbors */
		fNeighbors_1.RowAlias(facet, neighbors_1);
		fNeighbors_2.RowAlias(facet, neighbors_2);
		int nd1 = neighbors_1.Length();
		int nd2 = neighbors_2.Length();
	
		/* load functions */
		phi_1.Set(nip, nd1, fPhi_1(facet));
		phi_2.Set(nip, nd2, fPhi_2(facet));
			
		/* load derivatives */
		double* Dphi_1_ptr = fDPhi_1(facet);
		double* Dphi_2_ptr = fDPhi_2(facet);
		for (int i = 0; i < nip; i++)
		{
			Dphi_1[i].Set(nsd, nd1, Dphi_1_ptr);
			Dphi_2[i].Set(nsd, nd2, Dphi_2_ptr);

			Dphi_1_ptr += nsd*nd1;
			Dphi_2_ptr += nsd*nd2;
		}
	}
	else
	{
		cout << "\n MeshFreeSurfaceSupportT::LoadData: must use stored functions" << endl;
		throw ExceptionT::kGeneralFail;
	}
}

/* write MLS statistics */
void MeshFreeSurfaceSupportT::WriteStatistics(ostream& out) const
{
#pragma unused(out)
}

/*************************************************************************
* Private
*************************************************************************/

/* stores data for the given facet/side */
void MeshFreeSurfaceSupportT::StoreData(int facet, int side,
	ArrayT< AutoArrayT<int> >& neighbors, ArrayT< AutoArrayT<double> >& phi,
	ArrayT< AutoArrayT<double> >& Dphi)
{
	/* check */
	if (side != 0 && side != 1) throw ExceptionT::kOutOfRange;

	/* dimensions */
	int nip = neighbors.Length();	
	int nsd = fRefSurfaceShape.NumSD() + 1;

	/* destinations */
	AutoArrayT<int>& Count = (side == 0) ? fCount_1 : fCount_2;
	VariRaggedArray2DT<int>& Neighbors = (side == 0) ? fNeighbors_1 : fNeighbors_2;
	VariRaggedArray2DT<double>& Phi = (side == 0) ? fPhi_1 : fPhi_2;
	VariRaggedArray2DT<double>& DPhi = (side == 0) ? fDPhi_1 : fDPhi_2;

	/* check */
	if (facet < 0 || facet > Neighbors.MajorDim()) throw ExceptionT::kOutOfRange;

	/* collect all nodes used */
	fall_neighbors.Dimension(0);
	for (int i = 0; i < nip; i++)
		fall_neighbors.AppendUnique(neighbors[i]);
	int nnd = fall_neighbors.Length();

	/* generate local map */
	int shift;
	MakeInverseMap(fall_neighbors, fnode_map, shift);
		
	/* generate local data */
	int* map = fnode_map.Pointer();
	fphi.Dimension(nip*nnd);
	fDphi.Dimension(nip*nsd*nnd);
	
	fphi  = 0.0;
	fDphi = 0.0;
	for (int k = 0; k < nip; k++)
	{
		double*  p_phi_lhs = fphi.Pointer(k*nnd);
		double* p_Dphi_lhs = fDphi.Pointer(k*nsd*nnd);

		int nnd_k = neighbors[k].Length();
		int*    p_node = neighbors[k].Pointer();
		double*  p_phi = phi[k].Pointer();
		double* p_Dphi = Dphi[k].Pointer();
		for (int i = 0; i < nnd_k; i++)
		{
			int lnd = fnode_map[*p_node++ - shift];
		
			p_phi_lhs[lnd] = *p_phi++;

			double* lhs = p_Dphi_lhs + lnd;
			double* rhs = p_Dphi++;
			for (int j = 0; j < nsd; j++)
			{
				*lhs = *rhs;
				rhs += nnd_k;
				lhs += nnd;
			}
		}
	}

	/* write to database */
	if (facet == Neighbors.MajorDim())
	{
		Count.Append(fall_neighbors.Length());
		Neighbors.AddRow(facet, fall_neighbors);
		Phi.AddRow(facet, fphi);
		DPhi.AddRow(facet, fDphi);
	}
	else
	{
		Count[facet] = fall_neighbors.Length();
		Neighbors.SetRow(facet, fall_neighbors);
		Phi.SetRow(facet, fphi);
		DPhi.SetRow(facet, fDphi);
	}
}

/* make inverse map (filled with -1) */
void MeshFreeSurfaceSupportT::MakeInverseMap(const iAutoArrayT& map,
	AutoArrayT<int>& inv_map, int& shift) const
{
	/* range */
	int max;
	map.MinMax(shift, max);
	int range = max - shift + 1;
	
	/* dimension */
	inv_map.Dimension(range);
	inv_map = -1;

	/* make map */
	int dim = map.Length();
	for (int i = 0; i < dim; i++)
		inv_map[map[i] - shift] = i;	
}

/* returns the characteristic facet dimensions */
double MeshFreeSurfaceSupportT::FacetSize(const dArray2DT& facet_coords) const
{
	/* 3D surfaces */
	if (fRefSurfaceShape.NumSD() + 1 != 2)
	{
//TEMP
		const double* x1 = facet_coords(0);
		const double* x2 = facet_coords(1);
		double dx = x1[0] - x2[0];
		double dy = x1[1] - x2[1];
		double dz = x1[2] - x2[2];
		return sqrt(dx*dx + dy*dy + dz*dz);
	}
	/* 2D surfaces */
	else
	{
		const double* x1 = facet_coords(0);
		const double* x2 = facet_coords(1);
		double dx = x1[0] - x2[0];
		double dy = x1[1] - x2[1];
		return sqrt(dx*dx + dy*dy);
	}
}
