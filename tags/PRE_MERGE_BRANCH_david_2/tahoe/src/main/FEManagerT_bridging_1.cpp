/* $Id: FEManagerT_bridging_1.cpp,v 1.5 2004-07-22 08:41:38 paklein Exp $ */
#include "FEManagerT_bridging.h"
#ifdef BRIDGING_ELEMENT

#include "ModelManagerT.h"
#include "NodeManagerT.h"
#include "KBC_ControllerT.h"
#include "KBC_CardT.h"
#include "ofstreamT.h"
#include "ifstreamT.h"
#include "NLSolver.h"
#include "CommManagerT.h"

#include "BridgingScaleT.h"
#include "ParticleT.h"
#include "ParticlePairT.h"
#include "dSPMatrixT.h"
#include "EAMFCC3D.h"
#include "EAMT.h"
#include "ShapeFunctionT.h"
#include "dSymMatrixT.h"
#include "ElementSupportT.h"

/* headers needed to compute the correction for overlap */
#include "ContinuumElementT.h"
#include "SolidMatListT.h"
#include "FCC3D.h"
#include "Hex2D.h"
#include "Chain1D.h"
#include "BondLatticeT.h"
#include "nArrayGroupT.h"
#include "nVariMatrixT.h"

#include "LAdMatrixT.h" //TEMP
#include "CCSMatrixT.h"

/* debugging */
#define __DEBUG__ 1

/* atom/point types */
const char free_ = 'f';
const char not_free_ = 'n';

/* element types */
const char p_0 = 'a'; /* bond density = 0 */
const char p_1 = 'b'; /* bond density = 1 */
const char p_x = 'c'; /* unknown: 0 < bond density < 1 */

using namespace Tahoe;

/* compute internal correction for the overlap region */
void FEManagerT_bridging::CorrectOverlap_1(const RaggedArray2DT<int>& neighbors, const dArray2DT& coords, 
	double smoothing, double k2)
{
	const char caller[] = "FEManagerT_bridging::CorrectOverlap_1";

	/* collect nodes and cells in the overlap region */
	iArrayT overlap_cell;
	iArrayT overlap_node;
	CollectOverlapRegion_free(overlap_cell, overlap_node);
	if (overlap_node.Length() == 0) return;

	/* map of overlap_node in local numbering */
	InverseMapT overlap_node_map;
	overlap_node_map.SetOutOfRange(InverseMapT::MinusOne);
	overlap_node_map.SetMap(overlap_node);

	/* map of overlap_cell in local numbering */
	InverseMapT overlap_cell_map;
	overlap_cell_map.SetOutOfRange(InverseMapT::MinusOne);
	overlap_cell_map.SetMap(overlap_cell);

	/* compute reduced connectivity list */
	RaggedArray2DT<int> ghost_neighbors;
	GhostNodeBonds(neighbors, ghost_neighbors, overlap_cell_map);
	if (fPrintInput) {
		fMainOut << "\n Bonds to interpolation points (self as leading neighbor):\n";
		fMainOut << setw(kIntWidth) << "row" << "  n..." << '\n';
		iArrayT tmp(ghost_neighbors.Length(), ghost_neighbors.Pointer());
		tmp++;
		ghost_neighbors.WriteNumbered(fMainOut);
		tmp--;
		fMainOut.flush();
	}

	/* Cauchy-Born constitutive model */
	const ContinuumElementT* coarse = fFollowerCellData.ContinuumElement();
	const MaterialListT& mat_list = coarse->MaterialsList();
	if (mat_list.Length() > 1) 
		ExceptionT::GeneralFail(caller, "expecting only 1 material not %d", mat_list.Length());

	ContinuumMaterialT* cont_mat = mat_list[0];
	Hex2D* hex_2D = dynamic_cast<Hex2D*>(cont_mat);
	FCC3D* fcc_3D = dynamic_cast<FCC3D*>(cont_mat);
	Chain1D* mat_1D = dynamic_cast<Chain1D*>(cont_mat);
	if (!mat_1D && !hex_2D && !fcc_3D) ExceptionT::GeneralFail(caller, "could not resolve C-B material");

	const BondLatticeT& bond_lattice = (hex_2D) ? hex_2D->BondLattice() : ((fcc_3D) ? fcc_3D->BondLattice() : mat_1D->BondLattice());
	const dArray2DT& bonds = bond_lattice.Bonds();
	double V_0 = (hex_2D) ? hex_2D->CellVolume() : ((fcc_3D) ? fcc_3D->CellVolume() : mat_1D->CellVolume());
	double R_0 = (hex_2D) ? hex_2D->NearestNeighbor() : ((fcc_3D) ? fcc_3D->NearestNeighbor() : mat_1D->NearestNeighbor());

	/* unknown bond densities */
	dArray2DT bond_densities(overlap_cell.Length(), coarse->NumIP()*bonds.MajorDim());
	bond_densities = 1.0;

	/* works space */
	//TEMP
	//dArray2DT ddf_dpdp_i(overlap_cell.Length(), coarse->NumIP());
	LAdMatrixT ddf_dpdp_i;
	nVariMatrixT<double> ddf_dpdp_i_man(0, ddf_dpdp_i);
	dArray2DT p_i, dp_i, df_dp_i;
	nArray2DGroupT<double> ip_unknown_group(0, false, coarse->NumIP());
	ip_unknown_group.Register(p_i);
	ip_unknown_group.Register(dp_i);
	ip_unknown_group.Register(df_dp_i);

	/* solve unkowns */
	dArrayT sum_R_N(overlap_node.Length());
	dArrayT f_a(overlap_node.Length());
	dArrayT R_i(coords.MinorDim());
	AutoArrayT<int> overlap_cell_i;
	for (int i = 0; i < bonds.MajorDim(); i++) /* one bond density at a time */ {
		
		/* bond vector */
		R_i.SetToScaled(R_0, bonds(i));
		
		/* compute contribution from bonds to ghost atoms */
		ComputeSum_signR_Na(R_i, ghost_neighbors, coords, overlap_node_map, sum_R_N, overlap_cell_i);
		if (fPrintInput) /* debugging */ {
			
			/* coordinates */
			const NodeManagerT* node_manager = NodeManager();
			dArray2DT coords(overlap_node.Length(), node_manager->NumSD());
			coords.RowCollect(overlap_node, node_manager->InitialCoordinates());

			/* output file */
			StringT file;
			file.Root(fInputFile);
			file.Append(".RdN.b", i+1);
			file.Append(".out");
			
			/* write output */
			dArray2DT n_values(overlap_node.Length(), 1, sum_R_N.Pointer());
			ArrayT<StringT> n_labels(1);
			n_labels[0] = "f";
			WriteOutput(file, coords, overlap_node, n_values, n_labels);		
		}
		
		/* dimension the local problem */
		ddf_dpdp_i_man.SetDimensions(overlap_cell_i.Length()*coarse->NumIP());
		ip_unknown_group.SetMajorDimension(overlap_cell_i.Length(), false);

		/* initialize */
		p_i = 1.0;
		
		/* solve for bond densities */
		f_a = sum_R_N;
		Compute_df_dp_1(R_i, V_0, *coarse, overlap_cell_i, overlap_node_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
		double error_0 = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));
		cout << setw(kIntWidth) << i+1 << ": {" << R_i.no_wrap() << "}: " << error_0 << '\n';
		if (fPrintInput) {
			fMainOut << "overlap cells for bond: " << ": {" << R_i.no_wrap() << "}:\n";
			iArrayT tmp;
			tmp.Alias(overlap_cell_i);
			tmp++;
			fMainOut << tmp.wrap(5) << endl;
			tmp--;
		}
		double abs_tol = 1.0e-10;
		double rel_tol = 1.0e-10;
		double div_tol = 1.0e+06;		
		int max_iter = 1000;
		int iter = 0;
		double error = error_0;
		while (iter++ < max_iter && error > abs_tol && error/error_0 > rel_tol && error/error_0 < div_tol) {

			try {		
#if 0
			/* compute update vector (steepest descent) */
			for (int i = 0; i < p_i.Length(); i++) {
				if (fabs(ddf_dpdp_i(i,i)) < kSmall) 
					ExceptionT::BadJacobianDet(caller, "smoothing matrix is singular");
				dp_i[i] = -df_dp_i[i]/ddf_dpdp_i(i,i);
			}
#endif

			dp_i.SetToScaled(-1.0, df_dp_i);
			dArrayT tmp;
			tmp.Alias(dp_i);
#if __DEBUG__
			fMainOut << "f:\n" << tmp << endl;
#endif
			ddf_dpdp_i.LinearSolve(tmp);

			/* update densities */
			p_i += dp_i;
			
			/* recompute residual */			
			f_a = sum_R_N;
			Compute_df_dp_1(R_i, V_0, *coarse, overlap_cell_i, overlap_node_map, p_i, f_a, smoothing, k2, df_dp_i, ddf_dpdp_i);
			error = sqrt(dArrayT::Dot(df_dp_i,df_dp_i));
			cout << setw(kIntWidth) << iter << ": e/e_0 = " << error/error_0 << endl;
			
			} /* end try */
			catch (ExceptionT::CodeT error) {
				iter = max_iter; /* exit */
			}
		}
		
		if (error > abs_tol && error/error_0 > rel_tol) /* converged */ {
			cout << "FAIL" << endl;
			p_i = 1.0;
		}

		/* save result */
		int nb = bonds.MajorDim();
		for (int k = 0; k < overlap_cell_i.Length(); k++) {
			int cell = overlap_cell_map.Map(overlap_cell_i[k]);
			for (int j = 0; j < coarse->NumIP(); j++)
				bond_densities(cell, i+j*nb) = p_i(k,j);
		}
	}

	/* write densities */
	if (fPrintInput) {
	
		/* dimensions */
		const NodeManagerT* node_manager = NodeManager();
		int nsd = node_manager->NumSD();
		int nip = coarse->NumIP();
		int nen = coarse->NumElementNodes();

		/* shape functions */
		LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
		NodeManager()->RegisterCoordinates(element_coords);
		ShapeFunctionT shapes = ShapeFunctionT(coarse->ShapeFunction(), element_coords);
		shapes.Initialize();
	
		/* collect integration point coordinates */
		dArray2DT coords(overlap_cell.Length()*nip, nsd);
		dArrayT ip_coords;
		int index = 0;
		for (int i = 0; i < overlap_cell.Length(); i++) {
		
			/* collect element coordinates */
			const iArrayT& nodes = coarse->ElementCard(overlap_cell[i]).NodesX();
			element_coords.SetLocal(nodes);
			
			/* compute ip coordinates */
			for (int j = 0; j < nip; j++) {
				coords.RowAlias(index++, ip_coords);
				shapes.IPCoords(ip_coords, j);
			}
		}
		
		/* collect output labels */
		ArrayT<StringT> bond_labels(bonds.MajorDim());
		for (int i = 0; i < bond_labels.Length(); i++)
			bond_labels[i].Append("p_", i+1);
	
		/* output file */
		StringT file;
		file.Root(fInputFile);
		file.Append(".bond_density.out");
			
		/* write output */
		dArray2DT n_values(overlap_cell.Length()*nip, bonds.MajorDim(), bond_densities.Pointer()); /* maintain order, change shape */
		iArrayT point_map(coords.MajorDim());
		point_map.SetValueToPosition();
		WriteOutput(file, coords, point_map, n_values, bond_labels);			
	}

	/* write unknowns into the state variable space */
	ContinuumElementT* non_const_coarse = const_cast<ContinuumElementT*>(coarse);
	for (int i = 0; i < overlap_cell.Length(); i++) {
	
		/* element information */
		ElementCardT& element = non_const_coarse->ElementCard(overlap_cell[i]);
	
		/* allocate space */
		element.Dimension(0, bond_densities.MinorDim());
		
		/* copy in densities */
		element.DoubleData() = bond_densities(i);
	}
}

/* compute Cauchy-Born contribution to the nodal internal force */
void FEManagerT_bridging::Compute_df_dp_1(const dArrayT& R, double V_0, const ContinuumElementT& coarse, 
	const ArrayT<int>& overlap_cell, const InverseMapT& overlap_node_map, const dArray2DT& rho, 
	dArrayT& f_a, double smoothing, double k2, dArray2DT& df_dp, LAdMatrixT& ddf_dpdp) const
{
	/* dimensions */
	NodeManagerT* node_manager = NodeManager();
	int nsd = node_manager->NumSD();
	int nen = coarse.NumElementNodes();
	int nip = coarse.NumIP();

	/* element coordinates */
	LocalArrayT element_coords(LocalArrayT::kInitCoords, nen, nsd);
	node_manager->RegisterCoordinates(element_coords);

	/* shape functions */
	ShapeFunctionT shapes = ShapeFunctionT(coarse.ShapeFunction(), element_coords);
	shapes.Initialize();

	/* integrate bond density term over cells in overlap */
	dMatrixT grad_Na(nsd,nen);
	for (int i = 0; i < overlap_cell.Length(); i++) {
	
		/* set element information */
		const iArrayT& nodesX = coarse.ElementCard(overlap_cell[i]).NodesX();
		element_coords.SetLocal(nodesX);
		shapes.SetDerivatives();

		/* integration parameters */
		const double* p = rho(i);
		const double* j = shapes.IPDets();
		const double* w = shapes.IPWeights();
		
		/* integrate */
		shapes.TopIP();
		while (shapes.NextIP()) {

			/* get shape function gradients */
			shapes.GradNa(grad_Na);
		
			/* integration factor */
			double jw = (*j++)*(*w++);
			double pm1_jw_by_V = ((*p++) - 1.0)*jw/V_0;
		
			/* loop over nodes */
			for (int k = 0; k < nodesX.Length(); k++) {
			
				int index = overlap_node_map.Map(nodesX[k]);
				if (index > -1) /* node is in overlap */ {
				
					/* inner product of bond and shape function gradient */
					double R_dot_dN = grad_Na.DotCol(k, R);
				
					/* assemble */
					f_a[index] += (R_dot_dN*pm1_jw_by_V);
				}
			}
		}
	}

	/* gradient work space */
	const ParentDomainT& parent_domain = shapes.ParentDomain();	
	ArrayT<dMatrixT> ip_gradient(nip);
	dMatrixT jacobian_inv(nsd);
	dMatrixT A(nsd,nip), ATA(nip), ATA_int(nip);
	for (int i = 0; i < nip; i++) {
		ip_gradient[i].Dimension(nsd,nip);
		parent_domain.IPGradientTransform(i, ip_gradient[i]);
	}

	//TEMP - coupling all integration points
	dArray2DT df_a_dp(f_a.Length(), df_dp.Length());
	df_a_dp = 0.0;
	dArray2DT df_dp_smoothing = df_dp;
	df_dp_smoothing = 0.0;
#if 0
	dArrayT diag(nip);
	diag = 0.0;
#endif

	/* compute residual */
	dArray2DT df(nen, nip);
	//dArrayT dp_i_dp(nip);
	dMatrixT ddp_i_dpdp(nip);
	df_dp = 0.0;
	ddf_dpdp = 0.0;
	dArrayT element_rho;
	dArrayT element_force;
	for (int i = 0; i < overlap_cell.Length(); i++) {
	
		/* set element information */
		const iArrayT& nodesX = coarse.ElementCard(overlap_cell[i]).NodesX();
		element_coords.SetLocal(nodesX);
		shapes.SetDerivatives();

		/* integration parameters */
		const double* p = rho(i);
		const double* j = shapes.IPDets();
		const double* w = shapes.IPWeights();
		
		/* integrate */
		ATA_int = 0.0;
		ddp_i_dpdp = 0.0;
		df = 0.0;
		shapes.TopIP();
		while (shapes.NextIP()) {

			/* integration factor */
			int ip = shapes.CurrIP();
			double jw = (*j++)*(*w++);
			double pm1_jw = ((*p++) - 1.0)*jw;			
			double jw_by_V = jw/V_0;

			/* get shape function gradients */
			shapes.GradNa(grad_Na);

			/* integrate density gradient matrix */
			parent_domain.DomainJacobian(element_coords, ip, jacobian_inv);
			jacobian_inv.Inverse();
			A.MultATB(jacobian_inv, ip_gradient[ip]);
			ATA.MultATB(A,A);
			ATA_int.AddScaled(smoothing*jw, ATA);
			
			/* add penalized force term */
			df_dp(i,ip) += k2*pm1_jw;
			ddp_i_dpdp(ip,ip) += k2*jw;
					
			/* integrate the bond density term over the element */
			for (int k = 0; k < nodesX.Length(); k++) {
			
				int index = overlap_node_map.Map(nodesX[k]);
				if (index > -1) /* node is in overlap */ {
				
					/* inner product of bond and shape function gradient */
					double R_dot_dN = grad_Na.DotCol(k, R);
				
					/* assemble */
					df(k, ip) += (R_dot_dN*jw_by_V);
				}
			}
		}
		
		/* accumulate */
		for (int j = 0; j < nodesX.Length(); j++) {
			
			int index = overlap_node_map.Map(nodesX[j]);
			if (index > -1) /* node is in overlap */ {
			
				/* across all integration points */
				for (int k = 0; k < nip; k++) {
				
					/* force */
					df_dp(i,k) += f_a[index]*df(j,k);
					
					/* stiffness */
//					diag[k] += df(j,k)*df(j,k);
//					ddf_dpdp(i,k) += df(j,k)*df(j,k);
					
					//TEMP - accumulate
					df_a_dp(index, i*nip + k) += df(j,k);
				}
			}
		}
		
		/* regularization contribution to force and stiffness */
//		df_dp_smoothing.RowAlias(i, element_force);
		df_dp.RowAlias(i, element_force);
		rho.RowAlias(i, element_rho);		
		ATA_int.Multx(element_rho, element_force, 1.0, dMatrixT::kAccumulate);
//		for (int j = 0; j < nip; j++)
//			diag[j] += ATA_int(j,j);
//			ddf_dpdp(i,j) += ATA_int(j,j);

		/* penalty regularization */
		ATA_int += ddp_i_dpdp;
		
		ddf_dpdp.AddBlock(i*nip, i*nip, ATA_int);
	}

//cout << "df_dp: " << df_dp.no_wrap() << '\n';
//cout << "df_dp_smoothing: " << df_dp_smoothing.no_wrap() << '\n';
//df_dp += df_dp_smoothing;
//cout << "df_dp: " << df_dp.no_wrap() << '\n';

	//TEMP - add coupling term
	for (int i = 0; i < df_a_dp.MajorDim(); i++)
		ddf_dpdp.Outer(df_a_dp(i), df_a_dp(i), 1.0, dMatrixT::kAccumulate);

#if 0
	for (int i = 0; i < ddf_dpdp.Length(); i++)
		cout << ddf_dpdp[i] << '\n';
	cout.flush();
#endif
}

#endif  /* BRIDGING_ELEMENT */
