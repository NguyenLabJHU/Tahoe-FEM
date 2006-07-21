/* $Id: EAMT.cpp,v 1.68 2006-07-03 01:38:28 hspark Exp $ */

#include "EAMT.h"

#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"
#include "dSPMatrixT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"
#include "ParameterContainerT.h"

/* EAM potentials */
#include "ParadynEAMT.h"

using namespace Tahoe;

static int ipair = 1;
static int iEmb  = 1;

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */
const int kNumOutput = 8;
static const char* OutputNames[kNumOutput] = {
	"displacement",
	"potential_energy",
	"kinetic_energy",
	"stress",
	"strain",
	"slip_vector",	
	"centrosymmetry",
	"coordination_number"
};

/* constructor */
EAMT::EAMT(const ElementSupportT& support):
	ParticleT(support),
	fNeighbors(kMemoryHeadRoom),
	fNearestNeighbors(kMemoryHeadRoom),
	fRefNearestNeighbors(kMemoryHeadRoom),
	fEqnos(kMemoryHeadRoom),
	fForce_list_man(0, fForce_list),
	fElectronDensity_man(kMemoryHeadRoom, fElectronDensity, 1),
	fEmbeddingEnergy_man(kMemoryHeadRoom, fEmbeddingEnergy, 1),
	fEmbeddingForce_man(kMemoryHeadRoom, fEmbeddingForce, 1),
	fEmbeddingStiff_man(kMemoryHeadRoom, fEmbeddingStiff, 1),
	fElectronDensityMessageID(CommManagerT::kNULLMessageID),
	fEmbeddingEnergyMessageID(CommManagerT::kNULLMessageID),
	fEmbeddingForceMessageID(CommManagerT::kNULLMessageID),
	fEmbeddingStiffMessageID(CommManagerT::kNULLMessageID),
	frhop_rMessageID(CommManagerT::kNULLMessageID),
	fExternalEmbedForce(NULL),
	fExternalElecDensity(NULL),
	fExternalEmbedForceNodes(NULL),
	fExternalElecDensityNodes(NULL){
	SetName("particle_EAM");
}

/* collecting element group equation numbers */
void EAMT::Equations(AutoArrayT<const iArray2DT*>& eq_1,
		     AutoArrayT<const RaggedArray2DT<int>*>& eq_2)
{
#pragma unused(eq_1)

  /* dimension equations array */
  fEqnos.Configure(fNeighbors, NumDOF());

  /* get local equations numbers */
  Field().SetLocalEqnos(fNeighbors, fEqnos);

  /* add to list of equation numbers */
  eq_2.Append(&fEqnos);
}

/* collecting element geometry connectivities */
void EAMT::ConnectsX(AutoArrayT<const iArray2DT*>& connects) const
{
  /* NOTE: do not add anything to the geometry connectivity list */
#pragma unused(connects)
}

/* collecting element field connectivities */
void EAMT::ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
		     AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const
{
#pragma unused(connects_1)
  connects_2.AppendUnique(&fNeighbors);
}

void EAMT::WriteOutput(void)
{
  const char caller[] = "EAMT::WriteOutput";

  /* inherited */
  ParticleT::WriteOutput();
  
  /* muli-processor information */
  CommManagerT& comm_manager = ElementSupport().CommManager();
  const ArrayT<int>* proc_map = comm_manager.ProcessorMap();
  int rank = ElementSupport().Rank();

  /* dimensions */
  int ndof = NumDOF();
  int nstrs = dSymMatrixT::NumValues(ndof);

	/* number of output variables */
	iArrayT counts;
	SetOutputCount(fOutputFlags, counts);
	int num_output = counts.Sum();

	/* offsets to the different output values */
	iArrayT offsets(fOutputFlags.Length());
	offsets = 0;
	for (int i = 1; i < offsets.Length(); i++)
		offsets[i] = offsets[i-1] + counts[i-1];

  /* number of nodes */
  const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
  int non = (parition_nodes) ? parition_nodes->Length() : ElementSupport().NumNodes();

  /* map from partition node index */
  const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();

  /* output arrays length number of active nodes */
  dArray2DT n_values(non, num_output), e_values;
  n_values = 0.0;

  /* global coordinates */
  const dArray2DT& coords = ElementSupport().CurrentCoordinates();

  /* the field */
  const FieldT& field = Field();
  const dArray2DT& displacement = field[0];
  const dArray2DT* velocities = NULL;
  if (field.Order() > 0) velocities = &(field[1]);

  	/* atomic volume */
  	double V0 = 0.0;
  	if (NumSD() == 1)
    	V0 = fLatticeParameter;
	else if (NumSD() == 2)
		V0 = sqrt(3.0)*fLatticeParameter*fLatticeParameter/2.0; /* 2D hex */
	else /* 3D */
		V0 = fLatticeParameter*fLatticeParameter*fLatticeParameter/4.0; /* FCC */  
	
  /* collect mass per particle */
  int num_types = fTypeNames.Length();
  dArrayT mass(num_types);
  for (int i = 0; i < num_types; i++)
    mass[i] = fEAMProperties[fPropertiesMap(i,i)]->Mass();

	/* collect displacements */
  dArrayT vec, values_i;
  dSymMatrixT temp(ndof);
  for (int i = 0; i < non; i++) 
    {
      int   tag_i = (parition_nodes) ? (*parition_nodes)[i] : i;
      int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
      int  type_i = fType[tag_i];
      
      /* values for particle i */
      n_values.RowAlias(local_i, values_i);

		/* displacements */
		if (fOutputFlags[kDisplacement]) {
			vec.Set(ndof, values_i.Pointer(offsets[kDisplacement]));
			displacement.RowCopy(tag_i, vec);
		}

		/* kinetic contribution to the virial */
		if (velocities && fOutputFlags[kStress]) {
			velocities->RowAlias(tag_i, vec);
			temp.Outer(vec);
			int index = offsets[kStress];
		 	for (int cc = 0; cc < nstrs; cc++) {
		   		values_i[index++] = (fabs(V0) > kSmall) ? -mass[type_i]*temp[cc]/V0 : 0.0;
		 	}
		}
    }

	if (iEmb == 1 && (fOutputFlags[kPE] || fOutputFlags[kStress]))
    {
      /* get electron density */
      if (ndof == 2)
	  	GetRho2D(coords,fElectronDensity);
      else if (ndof == 3) 
	  	GetRho3D(coords,fElectronDensity);
	  else
	      ExceptionT::GeneralFail(caller);
	  
      /* exchange electron density information */
      comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);
      
      /* get embedding energy */
      GetEmbEnergy(coords,fElectronDensity,fEmbeddingEnergy);
	  
      /* exchange embedding energy information */
      comm_manager.AllGather(fEmbeddingEnergyMessageID, fEmbeddingEnergy);
    }

  /* EAM properties function pointers */
  EAMPropertyT::PairEnergyFunction  pair_energy_i = NULL;
  EAMPropertyT::PairEnergyFunction  pair_energy_j = NULL;

  EAMPropertyT::PairForceFunction  pair_force_i  = NULL;
  EAMPropertyT::PairForceFunction  pair_force_j  = NULL;

  EAMPropertyT::EDForceFunction ed_force_i = NULL;
  EAMPropertyT::EDForceFunction ed_force_j = NULL;
  fForce = 0.0;

	iArrayT neighbors;
	dArrayT x_i, x_j, r_ij(ndof);
  	dSymMatrixT vs_i(ndof);

  int current_property_i = -1;
  int current_property_j = -1;
		
	/* Loop i : run through neighbor list */
	for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* tags */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];		
      int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
      double* f_i = fForce(tag_i);

      /* values for particle i */
      n_values.RowAlias(local_i, values_i);		

		/* initialize the stress */
		if (fOutputFlags[kStress]) vs_i = 0.0;

		/* kinetic energy */
		if (velocities && fOutputFlags[kKE]) {
			velocities->RowAlias(tag_i, vec);
			values_i[offsets[kKE]] = 0.5*mass[type_i]*dArrayT::Dot(vec, vec);
		}

      coords.RowAlias(tag_i, x_i);

      /* Embedding Energy: E_i(rho_i) */
      if(iEmb == 1 && fOutputFlags[kPE]) values_i[offsets[kPE]] += fEmbeddingEnergy(tag_i,0);
		for (int j = 1; j < neighbors.Length(); j++)
		{

	  /* tags */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];		
	  double* f_j = fForce(tag_j);
			
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_energy_i  = fEAMProperties[property_i]->getPairEnergy();
	      pair_force_i  = fEAMProperties[property_i]->getPairForce();
	      ed_force_i    = fEAMProperties[property_i]->getElecDensForce();
	      current_property_i = property_i;
	    }
	 
	  int property_j = fPropertiesMap(type_j, type_i);
	  if (property_j != current_property_j)
	    {
	      pair_energy_j  = fEAMProperties[property_j]->getPairEnergy();
	      pair_force_j  = fEAMProperties[property_j]->getPairForce();
	      ed_force_j    = fEAMProperties[property_j]->getElecDensForce();
	      current_property_j = property_j;
	    }

	  /* global coordinates */
	  coords.RowAlias(tag_j, x_j);
	  
	  /* connecting vector */
	  r_ij.DiffOf(x_j, x_i);
	  double r = r_ij.Magnitude();

	  /* Pair Potential : phi = 0.5 * z_i z_j /r */
	  double phiby2 = 0.0;
	  if(ipair == 1 && fOutputFlags[kPE]) 
	  {
	    double z_i =  pair_energy_i(r, NULL, NULL);
	    double z_j =  pair_energy_j(r, NULL, NULL);
	    double phi =  z_i * z_j/r;
	    phiby2 = 0.5*phi;
	    values_i[offsets[kPE]] += phiby2;
	  }
	  
       	/* Compute Force  */
	  	double Fbyr=0.0;	  
		if (fOutputFlags[kStress]) {

			/* Component of force coming from Pair potential */
			if(ipair == 1) {
				double z_i = pair_energy_i(r,NULL,NULL);
				double z_j = pair_energy_j(r,NULL,NULL);
				double zp_i = pair_force_i(r,NULL,NULL);
				double zp_j = pair_force_j(r,NULL,NULL);
	      
				double E = z_i*z_j/r;
				double F = (z_i*zp_j + zp_i*z_j)/r - E/r;
	      
				Fbyr = F/r;
	    	}

			/* Component of force coming from Embedding energy */
			if(iEmb == 1) {
				double Ep_i   = fEmbeddingForce(tag_i,0);
				double Ep_j   = fEmbeddingForce(tag_j,0);
				double rhop_i = ed_force_i(r,NULL,NULL);
				double rhop_j = ed_force_j(r,NULL,NULL);
				
				double F =  Ep_j * rhop_i + Ep_i * rhop_j;
				Fbyr += F/r;
			}
	  		temp.Outer(r_ij);
			vs_i.AddScaled( 0.5*Fbyr,temp);
	  }

		/* second node may not be on processor */
		if (!proc_map || (*proc_map)[tag_j] == rank) 
		{
			int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;
			if (local_j < 0 || local_j >= n_values.MajorDim())
				cout << caller << ": out of range: " << local_j << '\n';
			else {

				/* potential energy */
				if (fOutputFlags[kPE]) n_values(local_j, offsets[kPE]) += phiby2;

				/* accumulate into stress into array */
				if (fOutputFlags[kStress]) {
					int index = offsets[kStress];
					for (int cc = 0; cc < nstrs; cc++)
						n_values(local_j, index++) += (fabs(V0) > kSmall) ? 0.5*Fbyr*temp[cc]/V0 : 0.0;
				}
			}
		}
		}

		/* copy stress into array */
		if (fOutputFlags[kStress]) {
			int index = offsets[kStress];
			for (int cc = 0; cc < nstrs; cc++)
				values_i[index++] += (fabs(V0) > kSmall) ? vs_i[cc]/V0 : 0.0;
		}
	}

	/* calculate centrosymmetry parameter */
	dArrayT csp;
	if (fOutputFlags[kCS]) {
		csp.Dimension(non);
		Calc_CSP(fNearestNeighbors, csp);
	}
	
	// calculate coordination number
	iArrayT cnarray;
	if (fOutputFlags[kCN]) {
		cnarray.Dimension(non);
		Calc_CN(fNearestNeighbors, cnarray);
	}

	/* slip vector, stress, and strain */
	dArray2DT s_values;
	if (fOutputFlags[kSlipVector] || fOutputFlags[kStress] || fOutputFlags[kStrain]) 
	{
		/* work space */
		int num_s_vals = nstrs + 1 + ndof; /* {strain, J, slip vector} */
		s_values.Dimension(non, num_s_vals);
		s_values = 0.0;

	    /* flag for specifying Lagrangian (0) or Eulerian (1) strain */ 
	    const int kEulerLagr = 0;

		/* calculate slip vector and strain */
		Calc_Slip_and_Strain(s_values, fRefNearestNeighbors, kEulerLagr);
	}

    /* combine strain, slip vector and centrosymmetry parameter into n_values list */
    for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
        /* row of neighbor list */
        fNeighbors.RowAlias(i, neighbors);

        /* tags */
        int   tag_i = neighbors[0]; /* self is 1st spot */
        int  type_i = fType[tag_i];
        int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

		/* strain */
		if (fOutputFlags[kStrain])
			for (int is = 0; is < nstrs; is++)
				n_values(local_i, offsets[kStrain]+is) = s_values(local_i,is);

		/* stress */
		if (fOutputFlags[kStress])
		{
			/* recover J, the determinant of the deformation gradient, for atom i
			 * and divide stress values by it */
			double J = s_values(local_i,nstrs);
	    	for (int is = 0; is < nstrs; is++)
	    		if (fabs(J) > kSmall)
	    			n_values(local_i, offsets[kStress]+is) /= J;
	    		else
	    			n_values(local_i, offsets[kStress]+is) = 0.0;
	    }

		/* slip vector */
		if (fOutputFlags[kSlipVector])
			for (int n = 0; n < ndof; n++)
				n_values(local_i, offsets[kSlipVector]+n) = s_values(local_i, nstrs+1+n);

		/* centrosymmetry */
		if (fOutputFlags[kCS])
			n_values(local_i, offsets[kCS]) = csp[local_i];
		
		// coordination number
		if (fOutputFlags[kCN])
			n_values(local_i, offsets[kCN]) = cnarray[local_i];
    }

	/* send */
	ElementSupport().WriteOutput(fOutputID, n_values, e_values);
}

/* compute the part of the stiffness matrix */
void EAMT::FormStiffness(const InverseMapT& col_to_col_eq_row_map,
			 const iArray2DT& col_eq, dSPMatrixT& stiffness)
{
  const char caller[] = "EAMT::FormStiffness";

  /* map should return -1 of out of range */
  if (col_to_col_eq_row_map.OutOfRange() != InverseMapT::MinusOne)
    ExceptionT::GeneralFail(caller, "inverse map out of range should return -1");

  /* assembly information */
  const ElementSupportT& support = ElementSupport();
  int group = Group();
  int ndof = NumDOF();
  fLHS.Dimension(2*ndof);
		
  /* global coordinates */
  const dArray2DT& coords = support.CurrentCoordinates();

  /* work space */
  dArrayT r_ij(NumDOF(), fRHS.Pointer());
  dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

  dArrayT r_ki(NumDOF()), r_kj(NumDOF());      

  /* run through neighbor list */
  const iArray2DT& field_eqnos = Field().Equations();
  iArrayT row_eqnos, col_eqnos; 
  iArrayT neighbors;
  dArrayT x_i, x_j, x_k;

  /* EAM properties function pointers */
  int current_property = -1;      
  EAMPropertyT::EDEnergyFunction ed_energy = NULL;

  EAMPropertyT::EDForceFunction ed_force_i = NULL;
  EAMPropertyT::EDForceFunction ed_force_j = NULL;

  fElectronDensity = 0.0;
  frhop_r = 0.0;
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* type */
      int  tag_i = neighbors[0]; /* self is 1st spot */
      int type_i = fType[tag_i];
      double* rp_i = frhop_r(tag_i);
		
      /* particle equations */
      field_eqnos.RowAlias(tag_i, row_eqnos);

      coords.RowAlias(tag_i, x_i);
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int tag_j = neighbors[j];
	  double* rp_j = frhop_r(tag_j);
			
	  /* particle is a target column */
	  int col_eq_index = col_to_col_eq_row_map.Map(tag_j);
	  if (col_eq_index != -1)
	    {
	      /* more particle info */
	      int type_j = fType[tag_j];

	      /* particle equations */
	      col_eq.RowAlias(col_eq_index, col_eqnos);

	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  ed_energy = fEAMProperties[property]->getElecDensEnergy();
		  current_property = property;
		}

	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);

	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();

	      fElectronDensity(tag_i,0) += ed_energy(r,NULL,NULL);	      
	      fElectronDensity(tag_j,0) += ed_energy(r,NULL,NULL);	      

	      double rhop_i = ed_force_i(r,NULL,NULL)/r;
	      double rhop_j = ed_force_j(r,NULL,NULL)/r;
	      for (int k = 0; k < ndof; k++)
		{
		  rp_i[k] +=  rhop_i * r_ij[k]; 
		  rp_j[k] += -rhop_j * r_ij[k]; 
		}
	    }
	}
    }

	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
 
  /* exchange electron density information */
  comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);

  /* exchange rhop * r information */
  comm_manager.AllGather(frhop_rMessageID, frhop_r);

  /* get embedding force */
  GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
  
  /* exchange embedding force information */
  comm_manager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

  /* get embedding stiffness */
  GetEmbStiff(coords,fElectronDensity,fEmbeddingStiff);
  
  /* exchange embedding stiffness information */
  comm_manager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);

  int current_property_i = -1;
  int current_property_j = -1;
  
  EAMPropertyT::PairEnergyFunction pair_energy_i = NULL;
  EAMPropertyT::PairEnergyFunction pair_energy_j = NULL;
  EAMPropertyT::PairForceFunction pair_force_i = NULL;
  EAMPropertyT::PairForceFunction pair_force_j = NULL;
  EAMPropertyT::PairStiffnessFunction pair_stiffness_i = NULL;
  EAMPropertyT::PairStiffnessFunction pair_stiffness_j = NULL;
  
  EAMPropertyT::EDStiffnessFunction ed_stiffness_i = NULL;
  EAMPropertyT::EDStiffnessFunction ed_stiffness_j = NULL;

  /* Loop i : run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* type */
      int  tag_i = neighbors[0]; /* self is 1st spot */
      int type_i = fType[tag_i];
      double* rp_i = frhop_r(tag_i);
		
      /* particle equations */
      field_eqnos.RowAlias(tag_i, row_eqnos);

      /* run though neighbors for one atom - first neighbor is self */
      coords.RowAlias(tag_i, x_i);

      /* Loop j */
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int tag_j = neighbors[j];
	  double* rp_j = frhop_r(tag_j);
			
	  /* particle is a target column */
	  int col_eq_index = col_to_col_eq_row_map.Map(tag_j);
	  if (col_eq_index != -1)
	    {
	      /* more particle info */
	      int type_j = fType[tag_j];

	      /* particle equations */
	      col_eq.RowAlias(col_eq_index, col_eqnos);

	      int property_i = fPropertiesMap(type_i, type_j);
	      if (property_i != current_property_i)
		{
		  pair_energy_i    = fEAMProperties[property_i]->getPairEnergy();
		  pair_force_i     = fEAMProperties[property_i]->getPairForce();
		  pair_stiffness_i = fEAMProperties[property_i]->getPairStiffness();

		  ed_force_i       = fEAMProperties[property_i]->getElecDensForce();
		  ed_stiffness_i   = fEAMProperties[property_i]->getElecDensStiffness();

		  current_property_i = property_i;
		}

	      int property_j = fPropertiesMap(type_j, type_i);
	      if (property_j != current_property_j)
		{
		  pair_energy_j    = fEAMProperties[property_j]->getPairEnergy();
		  pair_force_j     = fEAMProperties[property_j]->getPairForce();
		  pair_stiffness_j = fEAMProperties[property_j]->getPairStiffness();

		  ed_force_j       = fEAMProperties[property_j]->getElecDensForce();
		  ed_stiffness_j   = fEAMProperties[property_j]->getElecDensStiffness();

		  current_property_j = property_j;
		}


	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);


	      dArray2DT E_ij(ndof,ndof);
	      E_ij = 0.0;
	      if(iEmb == 1)
		{
		  /* Loop k */
		  for (int k = 1; k < neighbors.Length(); k++)
		    {
		      /* global tag */
		      int  tag_k = neighbors[k];
		      int type_k = fType[tag_k];

		      if(tag_k >  tag_j) 
			{
			  /* global coordinates */
			  coords.RowAlias(tag_k, x_k);
			  
			  r_ki.DiffOf(x_i, x_k);
			  r_kj.DiffOf(x_j, x_k);
			  
			  double rki = r_ki.Magnitude();
			  double rkj = r_kj.Magnitude();
			  
			  double EmbStiff = fEmbeddingStiff(tag_k,0) /rki /rkj ;
			  for (int m = 0; m < ndof; m++)
			    for (int n = 0; n < ndof; n++)
			      E_ij(m,n) += EmbStiff * r_ki[m] * r_kj[n]; 
			}
		    }
		}


	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
	      r_ji.SetToScaled(-1.0, r_ij);

	      /* Component of force coming from Pair potential */
	      if(ipair == 1)
		{
		  fLHS  = 0.0;
		  double z_i  = pair_energy_i(r, NULL, NULL);
		  double z_j  = pair_energy_j(r, NULL, NULL);
		  double zp_i = pair_force_i(r,NULL,NULL);
		  double zp_j = pair_force_j(r,NULL,NULL);
		  double zpp_i = pair_stiffness_i(r,NULL,NULL);
		  double zpp_j = pair_stiffness_j(r,NULL,NULL);
		  
		  double E =  z_i * z_j/r;
		  double F =  (z_i * zp_j + zp_i * z_j)/r - E/r;
		  double K =  (zpp_i*z_j + 2*zp_i * zp_j + z_i * zpp_j)/r - 2*F/r;
		  
		  double Fbyr = F/r;
		  
		  /* 1st term */
		  fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);
		
		  /* 2nd term */
		  fLHS.AddScaled(Fbyr, fOneOne);

		  /* assemble */
		  for (int p = 0; p < row_eqnos.Length(); p++)
		    for (int q = 0; q < col_eqnos.Length(); q++)
		      stiffness.AddElement(row_eqnos[p]-1, col_eqnos[q]-1, fLHS(p,q));

		}

	      /* Component of force coming from Embedding Energy */
	      if(iEmb == 1)
		{
		  fLHS  = 0.0;

		  double Ep_i   = fEmbeddingForce(tag_i,0); 
		  double Epp_i  = fEmbeddingStiff(tag_i,0); 
		  
		  double rhop_i  = ed_force_i(r,NULL,NULL);
		  double rhopp_i = ed_stiffness_i(r,NULL,NULL);
		  
		  double Ep_j   = fEmbeddingForce(tag_j,0); 
		  double Epp_j  = fEmbeddingStiff(tag_j,0); 
		  
		  double rhop_j  = ed_force_j(r,NULL,NULL);
		  double rhopp_j = ed_stiffness_j(r,NULL,NULL);
		  

		  double F =  Ep_j * rhop_i + Ep_i * rhop_j;
		  double Fbyr = F/r;

		  double K = Ep_j * rhopp_i + Ep_i * rhopp_j;

		  fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);
		  fLHS.AddScaled(Fbyr, fOneOne);

		  double T1 = Epp_i * rhop_j;
		  double T2 = Epp_j * rhop_i;
		  double L = rhop_i * rhop_j;

		  dArrayT El(2*ndof);
		  for (int k = 0; k < ndof; k++)
		    {
		      El[k]      =  (T1 * rp_i[k] + T2 * rp_j[k]);
		      El[k+ndof] = -(T1 * rp_i[k] + T2 * rp_j[k]);
		    }	
		  fLHS.Outer(fRHS, El, 1.0/r);
		  

		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < ndof; l++)
		      {
			fLHS(k,l)            += L * E_ij(k,l);
			fLHS(k+ndof,l)       += L * E_ij(k,l);
			fLHS(k,l+ndof)       += L * E_ij(k,l);
			fLHS(k+ndof,l+ndof)  += L * E_ij(k,l);
		      }

		  /* assemble */
		  for (int p = 0; p < row_eqnos.Length(); p++)
		    for (int q = 0; q < col_eqnos.Length(); q++)
		      stiffness.AddElement(row_eqnos[p]-1, col_eqnos[q]-1, fLHS(p,q));
		}
	    }
	}
    }
}

/* set external electron density pointers */
void EAMT::SetExternalElecDensity(const dArray2DT& elecdensity, const iArrayT& ghostatoms)
{
	fExternalElecDensity = &elecdensity;
	fExternalElecDensityNodes = &ghostatoms;
}

/* set external embedding force pointers */
void EAMT::SetExternalEmbedForce(const dArray2DT& embedforce, const iArrayT& ghostatoms)
{
	fExternalEmbedForce = &embedforce;
	fExternalEmbedForceNodes = &ghostatoms;
}

/* information about subordinate parameter lists */
void EAMT::DefineSubs(SubListT& sub_list) const
{
	/* inherited */
	ParticleT::DefineSubs(sub_list);

	/* interactions */
	sub_list.AddSub("EAM_particle_interaction", ParameterListT::OnePlus);

	/* output */
	sub_list.AddSub("particle_EAM_output", ParameterListT::ZeroOrOnce);
}

/* return the description of the given inline subordinate parameter list */
void EAMT::DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
	SubListT& sub_lists) const
{
	if (name == "EAM_property_choice")
	{
		order = ParameterListT::Choice;
		
		/* EAM potentials reading Paradyn parameters tables */
		sub_lists.AddSub("Paradyn_EAM");
	}
	else /* inherited */
		ParticleT::DefineInlineSub(name, order, sub_lists);
}

/* a pointer to the ParameterInterfaceT of the given subordinate */
ParameterInterfaceT* EAMT::NewSub(const StringT& name) const
{
	/* try to construct potential */
	EAMPropertyT* EAM_property = New_EAMProperty(name, false);
	if (EAM_property)
		return EAM_property;
	else if (name == "EAM_particle_interaction")
	{
		ParameterContainerT* interactions = new ParameterContainerT(name);
		interactions->SetSubSource(this);

		/* particle type labels */
		interactions->AddParameter(ParameterT::Word, "label_1");
		interactions->AddParameter(ParameterT::Word, "label_2");
	
		/* properties choice list */
		interactions->AddSub("EAM_property_choice", ParameterListT::Once, true);

		return interactions;
	}
	else if (name == "particle_EAM_output")
	{
		ParameterContainerT* output = new ParameterContainerT(name);
		
		/* all true by default */
		for (int i = 0; i < kNumOutput; i++) {
			ParameterT var(ParameterT::Integer, OutputNames[i]);
			var.SetDefault(1);
			output->AddParameter(var, ParameterListT::ZeroOrOnce);
		}

		return output;
	}		
	else /* inherited */
		return ParticleT::NewSub(name);
}

/* accept parameter list */
void EAMT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParticleT::TakeParameterList(list);

	/* output variables */
	fOutputFlags.Dimension(kNumOutput);
	fOutputFlags = 0;
	const ParameterListT* output = list.List("particle_EAM_output");
	if (output) 
	{
		/* set flags */
		for (int i = 0; i < kNumOutput; i++)
		{
			/* look for entry */
			const ParameterT* value = output->Parameter(OutputNames[i]);
			if (value) {
				int do_write = *value;
				if (do_write)
					fOutputFlags[i] = 1;
			}
		}
	}

	/* set the list of reference nearest neighbors */
	if (fOutputFlags[kSlipVector] || fOutputFlags[kStress] || fOutputFlags[kStrain])	
		SetRefNN(fNearestNeighbors, fRefNearestNeighbors);

	/* dimension */
	int ndof = NumDOF();
	fLHS.Dimension(2*ndof);
	fRHS.Dimension(2*ndof);

	/* constant matrix needed to calculate stiffness */
	fOneOne.Dimension(fLHS);
	dMatrixT one(ndof);
	one.Identity();
	fOneOne.SetBlock(0, 0, one);
	fOneOne.SetBlock(ndof, ndof, one);
	one *= -1;
	fOneOne.SetBlock(0, ndof, one);
	fOneOne.SetBlock(ndof, 0, one);
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* return number of values for each output variable */
void EAMT::SetOutputCount(const iArrayT& flags, iArrayT& counts) const
{
	/* dimension check */
	if (flags.Length() != kNumOutput)
		ExceptionT::SizeMismatch("EAMT::SetOutputCount");
	
	/* initialize */
	counts.Dimension(flags.Length());
	counts = 0;

	/* set output flags */
	if (flags[kDisplacement]) counts[kDisplacement] = NumDOF();
	if (flags[kPE]) counts[kPE] = 1;
	if (flags[kKE]) counts[kKE] = 1;
	if (flags[kCS]) counts[kCS] = 1;
	if (flags[kCN]) counts[kCN] = 1;
	if (flags[kStress]) counts[kStress] = dSymMatrixT::NumValues(NumSD());
	if (flags[kStrain]) counts[kStrain] = dSymMatrixT::NumValues(NumSD());
	if (flags[kSlipVector]) counts[kSlipVector] = NumDOF();
}

/* generate labels for output data */
void EAMT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
	const char caller[] = "EAMT::GenerateOutputLabels";
	int ndof = NumDOF();
	if (ndof > 3) ExceptionT::GeneralFail(caller);

	/* number of output variables */
	iArrayT counts;
	SetOutputCount(fOutputFlags, counts);
	int num_output = counts.Sum();

	/* offsets to the different output values */
	iArrayT offsets(fOutputFlags.Length());
	offsets = 0;
	for (int i = 1; i < offsets.Length(); i++)
		offsets[i] = offsets[i-1] + counts[i-1];

	/* initialize */
	labels.Dimension(num_output);

	/* displacement labels */
	if (fOutputFlags[kDisplacement]) {
		const char* disp[3] = {"D_X", "D_Y", "D_Z"};
		int index = offsets[kDisplacement];
		for (int i = 0; i < ndof; i++)
			labels[index++] = disp[i];
	}

	/* potential energy */
	if (fOutputFlags[kPE])
		labels[offsets[kPE]] = "PE";

	/* kinetic energy */
	if (fOutputFlags[kKE])
		labels[offsets[kKE]] = "KE";
	
	/* centrosymmetry */
	if (fOutputFlags[kCS])
		labels[offsets[kCS]] = "CS";
	
	// coordination number
	if (fOutputFlags[kCN])
		labels[offsets[kCN]] = "CN";
	
	/* slip vector */
	if (fOutputFlags[kSlipVector]) {
		const char* SV[3] = {"SV_X", "SV_Y", "SV_Z"};
		int index = offsets[kSlipVector];
		for (int i = 0; i < ndof; i++)
			labels[index++] = SV[i];
	}
	
	/* stress */
	if (fOutputFlags[kStress]) {
		const char* s1D[1] = {"s11"};
		const char* s2D[3] = {"s11", "s22", "s12"};
		const char* s3D[6] = {"s11", "s22", "s33", "s23", "s13", "s12"};

		const char** slabels = NULL;
		if (ndof == 1) slabels = s1D;
		else if (ndof == 2) slabels = s2D;
		else if (ndof == 3) slabels = s3D;
		else ExceptionT::GeneralFail(caller);	
		
		int nstrs = dSymMatrixT::NumValues(ndof);
		int index = offsets[kStress];
		for (int i = 0; i < nstrs; i++)
			labels[index++] = slabels[i];
	}
	
	/* strain */
	if (fOutputFlags[kStrain]) {		
		const char* e1D[1] = {"e11"};
		const char* e2D[3] = {"e11", "e12", "e22"};
		const char* e3D[6] = {"e11", "e12", "e13", "e22", "e23", "e33"};

		const char** elabels = NULL;
		if (ndof == 1) elabels = e1D;
		else if (ndof == 2) elabels = e2D;
		else if (ndof == 3) elabels = e3D;
		else ExceptionT::GeneralFail(caller);	
		
		int nstrn = dSymMatrixT::NumValues(ndof);
		int index = offsets[kStrain];
		for (int i = 0; i < nstrn; i++)
			labels[index++] = elabels[i];
	}
}

/* return a new EAM property or NULL if the name is invalid */
EAMPropertyT* EAMT::New_EAMProperty(const StringT& name, bool throw_on_fail) const
{
	if (name == "Paradyn_EAM")
		return new ParadynEAMT;
	else if (throw_on_fail) 
		ExceptionT::GeneralFail("EAMT::New_EAMProperty",
			"unrecognized potential \"%s\"", name.Pointer());

	return NULL;
}

/* form group contribution to the stiffness matrix */
void EAMT::LHSDriver(GlobalT::SystemTypeT sys_type)
{
  /* time integration parameters */
  double constK = 0.0;
  double constM = 0.0;
  int formK = fIntegrator->FormK(constK);
  int formM = fIntegrator->FormM(constM);

	/* assemble particle mass */
	if (formM) {
		
		/* collect mass per particle */
		int num_types = fTypeNames.Length();
		dArrayT mass(num_types);
		for (int i = 0; i < num_types; i++)
			mass[i] = fEAMProperties[fPropertiesMap(i,i)]->Mass();
		mass *= constM;

		AssembleParticleMass(mass);
	}

	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	
  /* assemble diagonal stiffness */
  if (formK && sys_type == GlobalT::kDiagonal)
    {
      /* assembly information */
      const ElementSupportT& support = ElementSupport();
      int group = Group();
      int ndof = NumDOF();
	
      /* global coordinates */
      const dArray2DT& coords = support.CurrentCoordinates();

      iArrayT neighbors;
      dArrayT x_i, x_j, r_ij(ndof);

      /* get electron density */
      if (ndof == 2) 
	  	GetRho2D(coords,fElectronDensity);
      else if (ndof == 3) 
	  	GetRho3D(coords,fElectronDensity);
	  else
	  	ExceptionT::GeneralFail();
	  
      /* exchange electron density information */
      comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);

      /* get embedding force */
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
	  
      /* exchange embedding force information */
      comm_manager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

      /* get embedding stiffness */
      GetEmbStiff(coords,fElectronDensity,fEmbeddingStiff);
      
	  /* exchange embedding stiffness information */
      comm_manager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);
   
      /* get rhop * r */
      frhop_r = 0.0;
      GetRhop_r(coords,frhop_r);
      /* exchange rhop * r information */
      comm_manager.AllGather(frhop_rMessageID, frhop_r);

      int current_property_i = -1;
      int current_property_j = -1;

      EAMPropertyT::PairEnergyFunction pair_energy_i = NULL;
      EAMPropertyT::PairEnergyFunction pair_energy_j = NULL;
      EAMPropertyT::PairForceFunction pair_force_i = NULL;
      EAMPropertyT::PairForceFunction pair_force_j = NULL;
      EAMPropertyT::PairStiffnessFunction pair_stiffness_i = NULL;
      EAMPropertyT::PairStiffnessFunction pair_stiffness_j = NULL;

      EAMPropertyT::EDForceFunction ed_force_i  = NULL;    
      EAMPropertyT::EDForceFunction ed_force_j  = NULL; 
      EAMPropertyT::EDStiffnessFunction ed_stiffness_i = NULL;
      EAMPropertyT::EDStiffnessFunction ed_stiffness_j = NULL;
      
      fForce = 0.0;
      /* Loop i : run through neighbor list */
      for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
	  /* row of neighbor list */
	  fNeighbors.RowAlias(i, neighbors);

	  /* type */
	  int  tag_i = neighbors[0]; /* self is 1st spot */
	  int type_i = fType[tag_i];
	  double* k_i = fForce(tag_i);
	  double* rp_i = frhop_r(tag_i);
		
	  coords.RowAlias(tag_i, x_i);

	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      double* k_j = fForce(tag_j);
	      double* rp_j = frhop_r(tag_j);
			
	      /* set EAM properties (if not already set) */
	      int property_i = fPropertiesMap(type_i, type_j);
	      if (property_i != current_property_i)
		{
		  pair_energy_i    = fEAMProperties[property_i]->getPairEnergy();
		  pair_force_i     = fEAMProperties[property_i]->getPairForce();
		  pair_stiffness_i = fEAMProperties[property_i]->getPairStiffness();

		  ed_force_i    = fEAMProperties[property_i]->getElecDensForce();
		  ed_stiffness_i= fEAMProperties[property_i]->getElecDensStiffness();

		  current_property_i = property_i;
		}	      

	      int property_j = fPropertiesMap(type_j, type_i);
	      if (property_j != current_property_j)
		{
		  pair_energy_j    = fEAMProperties[property_j]->getPairEnergy();
		  pair_force_j     = fEAMProperties[property_j]->getPairForce();
		  pair_stiffness_j = fEAMProperties[property_j]->getPairStiffness();

		  ed_force_j    = fEAMProperties[property_j]->getElecDensForce();
		  ed_stiffness_j= fEAMProperties[property_j]->getElecDensStiffness();

		  current_property_j = property_j;
		}
		
	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);
		
	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
			
	      /* Component of force coming from Pair potential */
	      if(ipair == 1)
		{
		  double z_i   = pair_energy_i(r, NULL, NULL);
		  double z_j   = pair_energy_j(r, NULL, NULL);
		  double zp_i  = pair_force_i(r,NULL,NULL);
		  double zp_j  = pair_force_j(r,NULL,NULL);
		  double zpp_i = pair_stiffness_i(r,NULL,NULL);
		  double zpp_j = pair_stiffness_j(r,NULL,NULL);
		  
		  double E = z_i * z_j/r;
		  double F = (z_i * zp_j + zp_i * z_j)/r - E/r;
		  double K = (zpp_i*z_j + 2*zp_i * zp_j + z_i * zpp_j)/r - 2*F/r;

		  double Fbyr = F/r;
		  
		  for (int k = 0; k < ndof; k++)
		    {
		      double r_k = r_ij[k]*r_ij[k]/r/r;
		      double K_k = constK*(K*r_k + Fbyr*(1.0 - r_k));
		      k_i[k] += K_k;
		      k_j[k] += K_k;
		    }
		}

	      /* Component of force coming from Embedding Energy */
	      if(iEmb == 1)
		{
		  double Ep_i   = fEmbeddingForce(tag_i,0); 
		  double Epp_i  = fEmbeddingStiff(tag_i,0); 

		  double rhop_i  = ed_force_i(r,NULL,NULL);
		  double rhopp_i = ed_stiffness_i(r,NULL,NULL);

		  double Ep_j   = fEmbeddingForce(tag_j,0); 
		  double Epp_j  = fEmbeddingStiff(tag_j,0); 

		  double rhop_j  = ed_force_j(r,NULL,NULL);
		  double rhopp_j = ed_stiffness_j(r,NULL,NULL);

		  double F =  Ep_j * rhop_i + Ep_i * rhop_j;
		  double Fbyr = F/r;

		  double K = Ep_j * rhopp_i + Ep_i * rhopp_j;

		  double T_i = Epp_j * rhop_i * rhop_i;
		  double T_j = Epp_i * rhop_j * rhop_j;

		  double L_i = Epp_i * rhop_j;
		  double L_j = Epp_j * rhop_i;

		  for (int k = 0; k < ndof; k++)
		    {
		      double r2_k = r_ij[k]/r;
		      double r_k  = r2_k * r2_k;
		      double K_k = K*r_k + Fbyr*(1.0 - r_k);
		      
		      k_i[k] += constK*K_k;
		      k_j[k] += constK*K_k;

		      double l_i = L_i * rp_i[k];
		      double l_j = L_j * rp_j[k];

		      k_i[k] += constK*(T_i * r_k + l_i * r2_k);
		      k_j[k] += constK*(T_j * r_k - l_j * r2_k);
		    }
		}
	    }
	}	

      /* assemble */
      support.AssembleLHS(group, fForce, Field().Equations());
    }
  else if (formK)
    {
      cout << "EAMT::LHSDriver, non-diagonal stiffness\n";
      /* assembly information */
      const ElementSupportT& support = ElementSupport();
      int group = Group();
      int ndof = NumDOF();
      fLHS.Dimension(2*ndof);

      /* global coordinates */
      const dArray2DT& coords = support.CurrentCoordinates();

      /* work space */
      dArrayT r_ij(NumDOF(), fRHS.Pointer());
      dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

      dArrayT r_ki(NumDOF()), r_kj(NumDOF());      

      const iArray2DT& field_eqnos = Field().Equations();
      iArray2DT pair_eqnos(2, ndof); 
      iArrayT pair(2);
      iArrayT neighbors;
      dArrayT x_i, x_j, x_k;

      /* EAM properties function pointers */
      if (ndof == 2) 
	  	GetRho2D(coords,fElectronDensity);
      else if(ndof == 3) 
	  	GetRho3D(coords,fElectronDensity);
	  else
	  	ExceptionT::GeneralFail();
		
      /* exchange electron density information */
      comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);

      /* get embedding force */
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
      
	  /* exchange embedding force information */
      comm_manager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

      /* get embedding stiffness */
      GetEmbStiff(coords,fElectronDensity,fEmbeddingStiff);
      
	  /* exchange embedding stiffness information */
      comm_manager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);
      
      /* get rhop * r */
      frhop_r = 0.0;
      GetRhop_r(coords,frhop_r);
      /* exchange rhop * r information */
      comm_manager.AllGather(frhop_rMessageID, frhop_r);

      int current_property_i = -1;
      int current_property_j = -1;

      EAMPropertyT::PairEnergyFunction pair_energy_i = NULL;
      EAMPropertyT::PairEnergyFunction pair_energy_j = NULL;
      EAMPropertyT::PairForceFunction pair_force_i = NULL;
      EAMPropertyT::PairForceFunction pair_force_j = NULL;
      EAMPropertyT::PairStiffnessFunction pair_stiffness_i = NULL;
      EAMPropertyT::PairStiffnessFunction pair_stiffness_j = NULL;

      EAMPropertyT::EDForceFunction ed_force_i = NULL;
      EAMPropertyT::EDForceFunction ed_force_j = NULL;
      EAMPropertyT::EDStiffnessFunction ed_stiffness_i = NULL;
      EAMPropertyT::EDStiffnessFunction ed_stiffness_j = NULL;

      /* Loop i : run through neighbor list */
      for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
	  /* row of neighbor list */
	  fNeighbors.RowAlias(i, neighbors);

	  /* type */
	  int  tag_i = neighbors[0]; /* self is 1st spot */
	  int type_i = fType[tag_i];
	  double* rp_i = frhop_r(tag_i);
	  pair[0] = tag_i;

	  coords.RowAlias(tag_i, x_i);

	  /* Loop j */
	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      double* rp_j = frhop_r(tag_j);

	      pair[1] = tag_j;
			
	      int property_i = fPropertiesMap(type_i, type_j);
	      if (property_i != current_property_i)
		{
		  pair_energy_i    = fEAMProperties[property_i]->getPairEnergy();
		  pair_force_i     = fEAMProperties[property_i]->getPairForce();
		  pair_stiffness_i = fEAMProperties[property_i]->getPairStiffness();

		  ed_force_i       = fEAMProperties[property_i]->getElecDensForce();
		  ed_stiffness_i   = fEAMProperties[property_i]->getElecDensStiffness();

		  current_property_i = property_i;
		}

	      int property_j = fPropertiesMap(type_j, type_i);
	      if (property_j != current_property_j)
		{
		  pair_energy_j    = fEAMProperties[property_j]->getPairEnergy();
		  pair_force_j     = fEAMProperties[property_j]->getPairForce();
		  pair_stiffness_j = fEAMProperties[property_j]->getPairStiffness();

		  ed_force_j       = fEAMProperties[property_j]->getElecDensForce();
		  ed_stiffness_j   = fEAMProperties[property_j]->getElecDensStiffness();

		  current_property_j = property_j;
		}
		
	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);

	      dArray2DT E_ij(ndof,ndof);
	      E_ij = 0.0;
	      if(iEmb == 1)
		{
		  /* Loop k */
		  for (int k = 1; k < neighbors.Length(); k++)
		    {
		      /* global tag */
		      int  tag_k = neighbors[k];
		      int type_k = fType[tag_k];

		      if(tag_k >  tag_j) 
			{
			  /* global coordinates */
			  coords.RowAlias(tag_k, x_k);
			  
			  r_ki.DiffOf(x_i, x_k);
			  r_kj.DiffOf(x_j, x_k);
			  
			  double rki = r_ki.Magnitude();
			  double rkj = r_kj.Magnitude();
			  
			  double EmbStiff = fEmbeddingStiff(tag_k,0) /rki /rkj ;
			  for (int m = 0; m < ndof; m++)
			    for (int n = 0; n < ndof; n++)
			      E_ij(m,n) += EmbStiff * r_ki[m] * r_kj[n]; 
			}
		    }
		}

	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
	      r_ji.SetToScaled(-1.0, r_ij);

	      /* Component of force coming from Pair potential */
	      if(ipair == 1)
	      {
		fLHS = 0.0;
		double z_i  = pair_energy_i(r, NULL, NULL);
		double z_j  = pair_energy_j(r, NULL, NULL);
		double zp_i = pair_force_i(r,NULL,NULL);
		double zp_j = pair_force_j(r,NULL,NULL);
		double zpp_i = pair_stiffness_i(r,NULL,NULL);
		double zpp_j = pair_stiffness_j(r,NULL,NULL);
		
		double E =  z_i * z_j/r;
		double F =  (z_i * zp_j + zp_i * z_j)/r - E/r;
		double K =  (zpp_i*z_j + 2*zp_i * zp_j + z_i * zpp_j)/r - 2*F/r;
		double Fbyr = F/r;

		/* 1st term */
		fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);

		/* 2nd term */
		fLHS.AddScaled(Fbyr, fOneOne);

		/* assemble */
		pair_eqnos.RowCollect(pair, field_eqnos);
		support.AssembleLHS(group, fLHS, pair_eqnos);
	      }
	      
	      /* Component of force coming from Embedding Energy */
	      if(iEmb == 1)
		{
		  fLHS  = 0.0;

		  double Ep_i   = fEmbeddingForce(tag_i,0); 
		  double Epp_i  = fEmbeddingStiff(tag_i,0); 
		  
		  double rhop_i  = ed_force_i(r,NULL,NULL);
		  double rhopp_i = ed_stiffness_i(r,NULL,NULL);
		  
		  double Ep_j   = fEmbeddingForce(tag_j,0); 
		  double Epp_j  = fEmbeddingStiff(tag_j,0); 
		  
		  double rhop_j  = ed_force_j(r,NULL,NULL);
		  double rhopp_j = ed_stiffness_j(r,NULL,NULL);
		  
		  double F =  Ep_j * rhop_i + Ep_i * rhop_j;
		  double Fbyr = F/r;

		  double K = Ep_j * rhopp_i + Ep_i * rhopp_j;

		  fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);
		  fLHS.AddScaled(Fbyr, fOneOne);

		  double T1 = Epp_i * rhop_j;
		  double T2 = Epp_j * rhop_i;
		  double L = rhop_i * rhop_j;

		  dArrayT El(2*ndof);
		  for (int k = 0; k < ndof; k++)
		    {
		      El[k]      =  (T1 * rp_i[k] + T2 * rp_j[k]);
		      El[k+ndof] = -(T1 * rp_i[k] + T2 * rp_j[k]);
		    }	
		  fLHS.Outer(fRHS, El, 1.0/r);
		  

		  for (int k = 0; k < ndof; k++)
		    for (int l = 0; l < ndof; l++)
		      {
			fLHS(k,l)            += L * E_ij(k,l);
			fLHS(k+ndof,l)       += L * E_ij(k,l);
			fLHS(k,l+ndof)       += L * E_ij(k,l);
			fLHS(k+ndof,l+ndof)  += L * E_ij(k,l);
		      }


		  /* assemble */
		  pair_eqnos.RowCollect(pair, field_eqnos);
		  support.AssembleLHS(group, fLHS, pair_eqnos);
		}
	    }
	} 
    }
}



/* form group contribution to the residual */
void EAMT::RHSDriver(void)
{
  int nsd = NumSD();
  if (nsd == 3)
    RHSDriver3D();
  else if (nsd == 2)
    RHSDriver2D();
  else
    ExceptionT::GeneralFail("EAMT::RHSDriver");
  
  ApplyDamping(fNeighbors);
	
  /* assemble */
  ElementSupport().AssembleRHS(Group(), fForce, Field().Equations());
}

void EAMT::RHSDriver2D(void)
{
  /* function name */
  const char caller[] = "EAMT::RHSDriver2D";

  /* check 2D */
  if (NumDOF() != 2) ExceptionT::GeneralFail(caller, "2D only: %d", NumDOF());

  /* time integration parameters */
  double constMa = 0.0;
  double constKd = 0.0;
  int formMa = fIntegrator->FormMa(constMa);
  int formKd = fIntegrator->FormKd(constKd);

  //TEMP - inertial force not implemented
  if (formMa) ExceptionT::GeneralFail(caller, "inertial force not implemented");

  /* assembly information */
  const ElementSupportT& support = ElementSupport();
  int group = Group();
  int ndof = NumDOF();
	
  /* global coordinates */
  const dArray2DT& coords = support.CurrentCoordinates();

	/* communication */
	CommManagerT& comm_manager = support.CommManager();

  if(iEmb == 1)
    {
      /* get electron density */
      GetRho2D(coords,fElectronDensity);
	  
      /* exchange electron density information */
      comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);
      
      /* get embedding force */
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);

      /* exchange embedding energy information */
      comm_manager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);
    }


  /* EAM properties function pointers */
  int current_property_i = -1;
  int current_property_j = -1;

  EAMPropertyT::PairEnergyFunction pair_energy_i = NULL;
  EAMPropertyT::PairEnergyFunction pair_energy_j = NULL;

  EAMPropertyT::PairForceFunction  pair_force_i  = NULL;
  EAMPropertyT::PairForceFunction  pair_force_j  = NULL;

  EAMPropertyT::EDForceFunction ed_force_i = NULL;
  EAMPropertyT::EDForceFunction ed_force_j = NULL;

  iArrayT neighbors;
  fForce = 0.0;

  /* Loop i: run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);
      
      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      double* f_i = fForce(tag_i);
      const double* x_i = coords(tag_i);

      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  double* f_j = fForce(tag_j);
	  const double* x_j = coords(tag_j);

	  /* set EAM property (if not already set) */
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_energy_i = fEAMProperties[property_i]->getPairEnergy();
	      pair_force_i  = fEAMProperties[property_i]->getPairForce();
	      ed_force_i    = fEAMProperties[property_i]->getElecDensForce();

	      current_property_i = property_i;
	    }

	  int property_j = fPropertiesMap(type_j, type_i);
	  if (property_j != current_property_j)
	    {
	      pair_energy_j = fEAMProperties[property_j]->getPairEnergy();
	      pair_force_j  = fEAMProperties[property_j]->getPairForce();

	      ed_force_j    = fEAMProperties[property_j]->getElecDensForce();
	      current_property_j = property_j;
	    }
			
	  /* Component of force coming from Pair potential */
	  if(ipair == 1)
	    {
	      double r_ij_0 = x_j[0] - x_i[0];
	      double r_ij_1 = x_j[1] - x_i[1];
	      double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);

	      double z_i = pair_energy_i(r,NULL,NULL);
	      double z_j = pair_energy_j(r,NULL,NULL);
	      double zp_i = pair_force_i(r,NULL,NULL);
	      double zp_j = pair_force_j(r,NULL,NULL);
	      
	      double E = z_i*z_j/r;
	      double F = (z_i*zp_j + zp_i*z_j)/r - E/r;
	      
	      double Fbyr = formKd*F/r;
	      
	      r_ij_0 *= Fbyr;
	      f_i[0] += r_ij_0;
	      f_j[0] +=-r_ij_0;
	      
	      r_ij_1 *= Fbyr;
	      f_i[1] += r_ij_1;
	      f_j[1] +=-r_ij_1;
	    }

	  /* Component of force coming from Embedding energy */
	  if(iEmb == 1)
	    {
	      double r_ij_0 = x_j[0] - x_i[0];
	      double r_ij_1 = x_j[1] - x_i[1];
	      double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);

	      double Ep_i   = fEmbeddingForce(tag_i,0);
	      double Ep_j   = fEmbeddingForce(tag_j,0);
	      double rhop_i = ed_force_i(r,NULL,NULL);
	      double rhop_j = ed_force_j(r,NULL,NULL);
	      
	      double F =  Ep_j * rhop_i + Ep_i * rhop_j;
	      double Fbyr = formKd*F/r;

	      r_ij_0 *= Fbyr;
	      f_i[0] +=  r_ij_0;
	      f_j[0] += -r_ij_0;
	      
	      r_ij_1 *= Fbyr;
	      f_i[1] +=  r_ij_1;
	      f_j[1] += -r_ij_1;
	    }
	}
    }
}

void EAMT::RHSDriver3D(void)
{
  /* function name */
  const char caller[] = "EAMT::RHSDriver3D";

  /* check 3D */
  if (NumDOF() != 3) ExceptionT::GeneralFail(caller, "3D only: %d", NumDOF());

  /* time integration parameters */
  double constMa = 0.0;
  double constKd = 0.0;
  int formMa = fIntegrator->FormMa(constMa);
  int formKd = fIntegrator->FormKd(constKd);

  //TEMP - inertial force not implemented
  if (formMa) ExceptionT::GeneralFail(caller, "inertial force not implemented");

  /* assembly information */
  const ElementSupportT& support = ElementSupport();
  int group = Group();
  int ndof = NumDOF();
	
  /* global coordinates */
  const dArray2DT& coords = support.CurrentCoordinates();

	/* communication */
	CommManagerT& comm_manager = support.CommManager();

  if(iEmb == 1)
    {
      /* get electron density */
      fElectronDensity = 0.0;
      GetRho3D(coords,fElectronDensity);
	  
	  if (fExternalElecDensity)
	  {
		dArrayT asdf(1);
		for (int i = 0; i < fExternalElecDensityNodes->Length(); i++)
		{
			fExternalElecDensity->RowAlias(i, asdf);
			fElectronDensity.SetRow((*fExternalElecDensityNodes)[i], asdf);
		}
	  }

      /* exchange electron density information */
      comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);

      /* get embedding force */
      fEmbeddingForce = 0.0;
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
	  
	  if (fExternalEmbedForce)
	  {
		dArrayT asdf(1);
		for (int i = 0; i < fExternalElecDensityNodes->Length(); i++)
		{
			fExternalEmbedForce->RowAlias(i, asdf);
			fEmbeddingForce.SetRow((*fExternalEmbedForceNodes)[i], asdf);
		}
	  }
	  
	  /* exchange embedding energy information */
      comm_manager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);
    }

 /* EAM properties function pointers */
  int current_property_i = -1;
  int current_property_j = -1;

  EAMPropertyT::PairEnergyFunction pair_energy_i = NULL;
  EAMPropertyT::PairEnergyFunction pair_energy_j = NULL;

  EAMPropertyT::PairForceFunction  pair_force_i  = NULL;
  EAMPropertyT::PairForceFunction  pair_force_j  = NULL;

  EAMPropertyT::EDForceFunction ed_force_i = NULL;
  EAMPropertyT::EDForceFunction ed_force_j = NULL;

  iArrayT neighbors;
  fForce = 0.0;

  /* Loop i : run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      double* f_i = fForce(tag_i);
      const double* x_i = coords(tag_i);

      /* Compute Force  */
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  double* f_j = fForce(tag_j);
	  const double* x_j = coords(tag_j);
		
	  /* set EAM property (if not already set) */
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_energy_i = fEAMProperties[property_i]->getPairEnergy();
	      pair_force_i  = fEAMProperties[property_i]->getPairForce();
	      ed_force_i    = fEAMProperties[property_i]->getElecDensForce();

	      current_property_i = property_i;
	    }
		
	  int property_j = fPropertiesMap(type_j, type_i);
	  if (property_j != current_property_j)
	    {
	      pair_energy_j = fEAMProperties[property_j]->getPairEnergy();
	      pair_force_j  = fEAMProperties[property_j]->getPairForce();
	      ed_force_j    = fEAMProperties[property_j]->getElecDensForce();

	      current_property_j = property_j;
	    }
		
	  /* Component of force coming from Pair potential */
	  if(ipair == 1)
	    {
	      double r_ij_0 = x_j[0] - x_i[0];
	      double r_ij_1 = x_j[1] - x_i[1];
	      double r_ij_2 = x_j[2] - x_i[2];
	      double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1 + r_ij_2*r_ij_2);

	      double z_i = pair_energy_i(r,NULL,NULL);
	      double z_j = pair_energy_j(r,NULL,NULL);
	      double zp_i = pair_force_i(r,NULL,NULL);
	      double zp_j = pair_force_j(r,NULL,NULL);
	      
	      double E = z_i*z_j/r;
	      double F = (z_i*zp_j + zp_i*z_j)/r - E/r;
	      double Fbyr = formKd*F/r;
		  
	      r_ij_0 *= Fbyr;
	      f_i[0] += r_ij_0;
	      f_j[0] +=-r_ij_0;
	      
	      r_ij_1 *= Fbyr;
	      f_i[1] += r_ij_1;
	      f_j[1] +=-r_ij_1;

	      r_ij_2 *= Fbyr;
	      f_i[2] += r_ij_2;
	      f_j[2] +=-r_ij_2;
	    }
	 	 

	  /* Component of force coming from Embedding energy */
	  if(iEmb == 1)
	    {
	      double r_ij_0 = x_j[0] - x_i[0];
	      double r_ij_1 = x_j[1] - x_i[1];
	      double r_ij_2 = x_j[2] - x_i[2];
	      double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1 + r_ij_2*r_ij_2);

	      double Ep_i   = fEmbeddingForce(tag_i,0);
	      double Ep_j   = fEmbeddingForce(tag_j,0);
	      double rhop_i = ed_force_i(r,NULL,NULL);
	      double rhop_j = ed_force_j(r,NULL,NULL);

	      double F =  Ep_j * rhop_i + Ep_i * rhop_j;
	      double Fbyr = formKd*F/r;

	      r_ij_0 *= Fbyr;
	      f_i[0] +=  r_ij_0;
	      f_j[0] += -r_ij_0;
	      
	      r_ij_1 *= Fbyr;
	      f_i[1] +=  r_ij_1;
	      f_j[1] += -r_ij_1;

	      r_ij_2 *= Fbyr;
	      f_i[2] +=  r_ij_2;
	      f_j[2] += -r_ij_2;	  
	 
	    }
	}
    }
}

/* set neighborlists */
void EAMT::SetConfiguration(void)
{
	/* inherited */
	ParticleT::SetConfiguration();

	/* reset neighbor lists */
	CommManagerT& comm_manager = ElementSupport().CommManager();
	const ArrayT<int>* part_nodes = comm_manager.PartitionNodes();
	if (fActiveParticles) 
		part_nodes = fActiveParticles;
	GenerateNeighborList(part_nodes, fNearestNeighborDistance, fNearestNeighbors, true, true);
	GenerateNeighborList(part_nodes, fNeighborDistance, fNeighbors, false, true);
	
	ofstreamT& out = ElementSupport().Output();
	out << "\n Neighbor statistics:\n";
	out << " Total number of neighbors . . . . . . . . . . . = "
	    << fNeighbors.Length() << '\n';
	out << " Minimum number of neighbors . . . . . . . . . . = " 
	    << fNeighbors.MinMinorDim(0) << '\n';
	out << " Maximum number of neighbors . . . . . . . . . . = " 
	    << fNeighbors.MaxMinorDim() << '\n';
	if (fNeighbors.MajorDim() > 0)
		out << " Average number of neighbors . . . . . . . . . . = " 
		    << double(fNeighbors.Length())/fNeighbors.MajorDim() << '\n';
	else
		out << " Average number of neighbors . . . . . . . . . . = " << 0 << '\n';

	/* verbose */
	if (ElementSupport().Logging() == GlobalT::kVerbose) {
		out << " Neighbor lists (self as leading neighbor):\n";
		out << setw(kIntWidth) << "row" << "  n..." << '\n';
		iArrayT tmp(fNeighbors.Length(), fNeighbors.Pointer());
		tmp++;
		fNeighbors.WriteNumbered(out);
		tmp--;
		out.flush();
	}

	/* initialize communications */
	if (fElectronDensityMessageID == CommManagerT::kNULLMessageID) {
		int ndof = NumDOF();
		fElectronDensityMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
		fEmbeddingEnergyMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
		fEmbeddingForceMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
		fEmbeddingStiffMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
		frhop_rMessageID = comm_manager.Init_AllGather(MessageT::Double, ndof);
		
		/* initialize memory manager */
		frhop_r_man.SetWard(kMemoryHeadRoom, frhop_r, NumDOF());
	}

	int nnd = ElementSupport().NumNodes();

	// ELECTRON DENSITY //
	/* reset the electron density array */
	fElectronDensity_man.SetMajorDimension(nnd, true);
  
	/* exchange type information */
	comm_manager.AllGather(fElectronDensityMessageID, fElectronDensity);

	// EMBEDDING ENERGY //
	/* reset the embedding energy array */
	fEmbeddingEnergy_man.SetMajorDimension(nnd, true);
  
	/* exchange type information */
	comm_manager.AllGather(fEmbeddingEnergyMessageID, fEmbeddingEnergy);

	/* reset the embedding force array */
	fEmbeddingForce_man.SetMajorDimension(nnd, true);
  
	/* exchange type information */
	comm_manager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

	/* reset the embedding stiffness array */
	fEmbeddingStiff_man.SetMajorDimension(nnd, true);
  
	/* exchange type information */
	comm_manager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);

	// OTHER  //
	/* reset the rhop * r array */
	frhop_r_man.SetMajorDimension(nnd, true);
  
	/* exchange type information */
	comm_manager.AllGather(frhop_rMessageID, frhop_r);
}

/* extract the properties information from the parameter list. See ParticleT::ExtractProperties */
void EAMT::ExtractProperties(const ParameterListT& list, const ArrayT<StringT>& type_names,
	ArrayT<ParticlePropertyT*>& properties, nMatrixT<int>& properties_map)
{
	const char caller[] = "EAMT::ExtractProperties";

	/* check number of interactions */
	int num_props = list.NumLists("EAM_particle_interaction");
	int dim = 0;
	for (int i = 0; i < properties_map.Rows(); i++)
		dim += properties_map.Rows() - i;
	if (dim != num_props)
		ExceptionT::GeneralFail(caller, "%d types requires %d \"EAM_particle_interaction\"",
			properties_map.Rows(), dim);

	/* read properties */
	fEAMProperties.Dimension(num_props);
	for (int i = 0; i < num_props; i++) {

		const ParameterListT& interaction = list.GetList("EAM_particle_interaction", i);
		
		/* type names */
		const StringT& label_1 = interaction.GetParameter("label_1");
		const StringT& label_2 = interaction.GetParameter("label_2");
		
		/* resolve index */
		int index_1 = -1, index_2 = -1;
		for (int j = 0; index_1 == -1 && j < type_names.Length(); j++)
			if (type_names[j] == label_1) index_1 = j;
		if (index_1 == -1) ExceptionT::GeneralFail(caller, "could not resolve index of \"%s\"", label_1.Pointer());

		for (int j = 0; index_2 == -1 && j < type_names.Length(); j++)
			if (type_names[j] == label_2) index_2 = j;
		if (index_2 == -1) ExceptionT::GeneralFail(caller, "could not resolve index of \"%s\"", label_2.Pointer());
		
		/* set properies map */
		if (properties_map(index_1, index_2) != -1)
			ExceptionT::GeneralFail(caller, "%s-%s interaction is already defined",
				label_1.Pointer(), label_2.Pointer());
		properties_map(index_1, index_2) = properties_map(index_2, index_1) = i; /* symmetric */
		
		/* read property */
		const ParameterListT& property = interaction.GetListChoice(*this, "EAM_property_choice");
		EAMPropertyT* EAM_prop = New_EAMProperty(property.Name(), true);
		EAM_prop->TakeParameterList(property);
		fEAMProperties[i] = EAM_prop;
	}
	
	/* copy */
	properties.Dimension(fEAMProperties.Length());
	for (int i = 0; i < properties.Length(); i++)
		properties[i] = fEAMProperties[i];
}

void EAMT::GetRho2D(const dArray2DT& coords,dArray2DT& rho)
{
  int current_property = -1;
  EAMPropertyT::EDEnergyFunction ed_energy = NULL;
  iArrayT neighbors;
  rho = 0.0;

  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);
      
      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      const double* x_i = coords(tag_i);

      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  const double* x_j = coords(tag_j);

	  int property = fPropertiesMap(type_i, type_j);
	  if (property != current_property)
	    {
	      ed_energy  = fEAMProperties[property]->getElecDensEnergy();
	      current_property = property;
	    }

	  /* connecting vector */
	  double r_ij_0 = x_j[0] - x_i[0];
	  double r_ij_1 = x_j[1] - x_i[1];
	  double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);
		
	  rho(tag_i,0) += ed_energy(r,NULL,NULL); 
	  rho(tag_j,0) += ed_energy(r,NULL,NULL);
	}
    }
}


void EAMT::GetRho3D(const dArray2DT& coords,dArray2DT& rho)
{
  int current_property = -1;
  EAMPropertyT::EDEnergyFunction ed_energy = NULL;
  iArrayT neighbors;
  rho = 0.0;

  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);
      
      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      const double* x_i = coords(tag_i);

      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  const double* x_j = coords(tag_j);

	  int property = fPropertiesMap(type_i, type_j);
	  if (property != current_property)
	    {
	      ed_energy  = fEAMProperties[property]->getElecDensEnergy();
	      current_property = property;
	    }
		
	  /* connecting vector */
	  double r_ij_0 = x_j[0] - x_i[0];
	  double r_ij_1 = x_j[1] - x_i[1];
	  double r_ij_2 = x_j[2] - x_i[2];
	  double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1 + r_ij_2*r_ij_2);
	  
	  rho(tag_i,0) += ed_energy(r,NULL,NULL); 
	  rho(tag_j,0) += ed_energy(r,NULL,NULL); 
	}
    }
}


void EAMT::GetRhop_r(const dArray2DT& coords,dArray2DT& rho)
{

  int ndof = NumDOF();
  int current_property_i = -1;
  int current_property_j = -1;
  
  EAMPropertyT::EDForceFunction ed_force_i  = NULL;    
  EAMPropertyT::EDForceFunction ed_force_j  = NULL; 

  iArrayT neighbors;
  dArrayT x_i, x_j, r_ij(ndof);
  
  rho = 0.0;
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      fNeighbors.RowAlias(i, neighbors);

	  int  tag_i = neighbors[0]; 
	  int type_i = fType[tag_i];
	  double* rp_i = rho(tag_i);
		
	  coords.RowAlias(tag_i, x_i);

	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      double* rp_j = rho(tag_j);
	      
	      coords.RowAlias(tag_j, x_j);
			
	      int property_i = fPropertiesMap(type_i, type_j);
	      if (property_i != current_property_i)
		{
		  ed_force_i  = fEAMProperties[property_i]->getElecDensForce();
		  current_property_i = property_i;
		}

	      int property_j = fPropertiesMap(type_j, type_i);
	      if (property_j != current_property_j)
		{
		  ed_force_j  = fEAMProperties[property_j]->getElecDensForce();
		  current_property_j = property_j;
		}	      	      
		
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
	      double rhop_i = ed_force_i(r,NULL,NULL)/r;
	      double rhop_j = ed_force_j(r,NULL,NULL)/r;
	      for (int k = 0; k < ndof; k++)
		{
		  rp_i[k] +=  rhop_j * r_ij[k]; 
		  rp_j[k] += -rhop_i * r_ij[k]; 
		}
	    }
	}
}


void EAMT::GetEmbEnergy(const dArray2DT& coords,const dArray2DT rho,
		  dArray2DT& Emb)
{
#pragma unused(coords)

  int current_property = -1;
  EAMPropertyT::EmbedEnergyFunction emb_energy = NULL;
  iArrayT neighbors;
  Emb = 0.0;
  
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      fNeighbors.RowAlias(i, neighbors);
      
      int   tag_i = neighbors[0]; 
      int  type_i = fType[tag_i];
      
      int property = fPropertiesMap(type_i, type_i);
      if (property != current_property)
	{
	  emb_energy  = fEAMProperties[property]->getEmbedEnergy();
	  current_property = property;
	}

      Emb(tag_i,0) = emb_energy(rho(tag_i,0),NULL,NULL); 
    }
}


void EAMT::GetEmbForce(const dArray2DT& coords,const dArray2DT rho,
	  	       dArray2DT& Emb)
{
#pragma unused(coords)

  int current_property = -1;
  EAMPropertyT::EmbedForceFunction emb_force = NULL;
  iArrayT neighbors;
  Emb = 0.0;
  
    for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      fNeighbors.RowAlias(i, neighbors);
      
      int   tag_i = neighbors[0]; 
      int  type_i = fType[tag_i];
	  
	  int property = fPropertiesMap(type_i, type_i);
      if (property != current_property)
	  {
		emb_force  = fEAMProperties[property]->getEmbedForce();
		current_property = property;
	  }
      //emb_force  = fEAMProperties[type_i]->getEmbedForce();
      Emb(tag_i,0) = emb_force(rho(tag_i,0),NULL,NULL); 
    }  
}

void EAMT::GetEmbStiff(const dArray2DT& coords,const dArray2DT rho,
		       dArray2DT& Emb)
{
#pragma unused(coords)

  int current_property = -1;
  EAMPropertyT::EmbedStiffnessFunction emb_stiffness = NULL;
  iArrayT neighbors;
  Emb = 0.0;
  
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      fNeighbors.RowAlias(i, neighbors);
      
      int   tag_i = neighbors[0]; 
      int  type_i = fType[tag_i];
      
      int property = fPropertiesMap(type_i, type_i);
      if (property != current_property)
	{
	  emb_stiffness  = fEAMProperties[property]->getEmbedStiffness();
	  current_property = property;
	}
 
      Emb(tag_i,0) = emb_stiffness(rho(tag_i,0),NULL,NULL); 
    }  
}
