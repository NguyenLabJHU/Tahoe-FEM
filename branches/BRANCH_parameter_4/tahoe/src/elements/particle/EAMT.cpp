/* $Id: EAMT.cpp,v 1.61 2004-06-28 22:41:39 hspark Exp $ */

#include "EAMT.h"

#include "ifstreamT.h"
#include "ofstreamT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"
#include "dSPMatrixT.h"
#include "dSymMatrixT.h"
#include "dArray2DT.h"

/* EAM potentials */
#include "ParadynEAMT.h"

using namespace Tahoe;

static int ipair = 1;
static int iEmb  = 1;

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */

/* constructor */
EAMT::EAMT(const ElementSupportT& support, const FieldT& field):
  ParticleT(support, field),
  fNeighbors(kMemoryHeadRoom),
  NearestNeighbors(kMemoryHeadRoom),
  RefNearestNeighbors(kMemoryHeadRoom),
  fEqnos(kMemoryHeadRoom),
  fForce_list_man(0, fForce_list),
  fElectronDensity_man(kMemoryHeadRoom, fElectronDensity, 1),
  fEmbeddingEnergy_man(kMemoryHeadRoom, fEmbeddingEnergy, 1),
  fEmbeddingForce_man(kMemoryHeadRoom, fEmbeddingForce, 1),
  fEmbeddingStiff_man(kMemoryHeadRoom, fEmbeddingStiff, 1),
  frhop_r_man(kMemoryHeadRoom, frhop_r,NumDOF()),
  fExternalEmbedForce(NULL),
  fExternalElecDensity(NULL),
  fExternalEmbedForceNodes(NULL),
  fExternalElecDensityNodes(NULL)
{
	SetName("particle_eam");
}

EAMT::EAMT(const ElementSupportT& support):
  ParticleT(support),
  fNeighbors(kMemoryHeadRoom),
  NearestNeighbors(kMemoryHeadRoom),
  RefNearestNeighbors(kMemoryHeadRoom),
  fEqnos(kMemoryHeadRoom),
  fForce_list_man(0, fForce_list),
  fElectronDensity_man(kMemoryHeadRoom, fElectronDensity, 1),
  fEmbeddingEnergy_man(kMemoryHeadRoom, fEmbeddingEnergy, 1),
  fEmbeddingForce_man(kMemoryHeadRoom, fEmbeddingForce, 1),
  fEmbeddingStiff_man(kMemoryHeadRoom, fEmbeddingStiff, 1),
  fExternalEmbedForce(NULL),
  fExternalElecDensity(NULL),
  fExternalEmbedForceNodes(NULL),
  fExternalElecDensityNodes(NULL)
{
	SetName("particle_eam");
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

/* class initialization */
void EAMT::Initialize(void)
{
  cout << "Initialization Phase\n";
  
	/* muli-processor information */
	CommManagerT& comm_manager = ElementSupport().CommManager();

  /* set up communication of electron density information */
  fElectronDensityMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
  fElectronDensity_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of embedding energy information */
  fEmbeddingEnergyMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
  fEmbeddingEnergy_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of embedding force information */
  fEmbeddingForceMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
  fEmbeddingForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of embedding stiffness information */
  fEmbeddingStiffMessageID = comm_manager.Init_AllGather(MessageT::Double, 1);
  fEmbeddingStiff_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of rhop * r information */
  frhop_rMessageID = comm_manager.Init_AllGather(MessageT::Double, NumDOF());
  frhop_r_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* inherited */
  ParticleT::Initialize();

  ParticleT::SetRefNN(NearestNeighbors,RefNearestNeighbors);

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
  int num_output = ndof + 2; /* displacement + PE + KE */

#ifndef NO_PARTICLE_STRESS_OUTPUT
  num_output++; /*includes centrosymmetry*/
  num_output+=ndof; /*some more for slip vector*/
#endif /* NO_PARTICLE_STRESS_OUTPUT */

  /* number of nodes */
  const ArrayT<int>* parition_nodes = comm_manager.PartitionNodes();
  int non = (parition_nodes) ? parition_nodes->Length() : ElementSupport().NumNodes();

  /* map from partition node index */
  const InverseMapT* inverse_map = comm_manager.PartitionNodes_inv();

#ifndef NO_PARTICLE_STRESS_OUTPUT
  dSymMatrixT vs_i(ndof), temp(ndof);
  int num_stresses=vs_i.NumValues(ndof);
  //dArray2DT vsvalues(non, num_stresses);
  num_output +=num_stresses; 
  num_output += num_stresses; //another for the strain
#endif

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
  dArrayT mass(fNumTypes);
  for (int i = 0; i < fNumTypes; i++)
    mass[i] = fEAMProperties[fPropertiesMap(i,i)]->Mass();

 /* collect displacements */
  dArrayT vec, values_i;
  for (int i = 0; i < non; i++) 
    {
      int   tag_i = (parition_nodes) ? (*parition_nodes)[i] : i;
      int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
      int  type_i = fType[tag_i];
      
      /* values for particle i */
      n_values.RowAlias(local_i, values_i);
      
      /* copy in */
      vec.Set(ndof, values_i.Pointer());
      displacement.RowCopy(tag_i, vec);
  
#ifndef NO_PARTICLE_STRESS_OUTPUT
	/* kinetic contribution to the virial */
		if (velocities) {
			velocities->RowAlias(tag_i, vec);
			temp.Outer(vec);
		 	for (int cc = 0; cc < num_stresses; cc++) {
				int ndex = ndof+2+cc;
		   		values_i[ndex] = (fabs(V0) > kSmall) ? -mass[type_i]*temp[cc]/V0 : 0.0;
		 	}
		}
#endif /* NO_PARTICLE_STRESS_OUTPUT */
    }
  if(iEmb == 1)
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

  int current_property_i = -1;
  int current_property_j = -1;
		
  /* Loop i : run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

#ifndef NO_PARTICLE_STRESS_OUTPUT
      vs_i=0.0;
#endif /* NO_PARTICLE_STRESS_OUTPUT */

      /* tags */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];		
      int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;
      double* f_i = fForce(tag_i);

      /* values for particle i */
      n_values.RowAlias(local_i, values_i);		

      /* kinetic energy */
      if (velocities)
	{
	  velocities->RowAlias(tag_i, vec);
	  values_i[ndof+1] = 0.5*mass[type_i]*dArrayT::Dot(vec, vec);
		    
	}
      coords.RowAlias(tag_i, x_i);

      /* Embedding Energy: E_i(rho_i) */
      if(iEmb == 1) values_i[ndof] += fEmbeddingEnergy(tag_i,0);
      	  
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
	  if(ipair == 1) 
	  {
	    double z_i =  pair_energy_i(r, NULL, NULL);
	    double z_j =  pair_energy_j(r, NULL, NULL);
	    double phi =  z_i * z_j/r;
	    phiby2 = 0.5*phi;
	  }
	  values_i[ndof] +=  phiby2;
	  
       /* Compute Force  */
	  double Fbyr=0.0;
	  /* Component of force coming from Pair potential */
	  if(ipair == 1)
	    {
	      double z_i = pair_energy_i(r,NULL,NULL);
	      double z_j = pair_energy_j(r,NULL,NULL);
	      double zp_i = pair_force_i(r,NULL,NULL);
	      double zp_j = pair_force_j(r,NULL,NULL);
	      
	      double E = z_i*z_j/r;
	      double F = (z_i*zp_j + zp_i*z_j)/r - E/r;
	      
	      Fbyr = F/r;
	    }

	  /* Component of force coming from Embedding energy */
	      if(iEmb == 1){
	    
	      double Ep_i   = fEmbeddingForce(tag_i,0);
	      double Ep_j   = fEmbeddingForce(tag_j,0);
	      double rhop_i = ed_force_i(r,NULL,NULL);
	      double rhop_j = ed_force_j(r,NULL,NULL);

	      double F =  Ep_j * rhop_i + Ep_i * rhop_j;
	      Fbyr += F/r;
	      }

#ifndef NO_PARTICLE_STRESS_OUTPUT
		temp.Outer(r_ij);
		vs_i.AddScaled( 0.5*Fbyr,temp);
#endif /* NO_PARTICLE_STRESS_OUTPUT */

		/* second node may not be on processor */
		if (!proc_map || (*proc_map)[tag_j] == rank) 
		{
			int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;
			if (local_j < 0 || local_j >= n_values.MajorDim())
				cout << caller << ": out of range: " << local_j << '\n';
			else {

				/* potential energy */
				n_values(local_j, ndof) += phiby2;

#ifndef NO_PARTICLE_STRESS_OUTPUT
				/* accumulate into stress into array */
				for (int cc = 0; cc < num_stresses; cc++) {
					int ndex = ndof+2+cc;
					n_values(local_j, ndex) += (fabs(V0) > kSmall) ? 0.5*Fbyr*temp[cc]/V0 : 0.0;
				}
#endif /* NO_PARTICLE_STRESS_OUTPUT */
			}	  
		}
	}

#ifndef NO_PARTICLE_STRESS_OUTPUT
	          /* copy stress into array */
	          for (int cc = 0; cc < num_stresses; cc++) {
	            int ndex = ndof+2+cc;
                values_i[ndex] += (fabs(V0) > kSmall) ? vs_i[cc]/V0 : 0.0;
	          }
#endif
	}

#ifndef NO_PARTICLE_STRESS_OUTPUT
    int num_s_vals = num_stresses+1+ndof+1;
	dArray2DT s_values(non,num_s_vals);

    /* flag for specifying Lagrangian (0) or Eulerian (1) strain */
    const int kEulerLagr = 0;
    /* calculate slip vector and strain */
    Calc_Slip_and_Strain(s_values,RefNearestNeighbors,kEulerLagr);
    /* calculate centrosymmetry parameter */
    Calc_CSP(s_values, NearestNeighbors);

    /* combine strain, slip vector and centrosymmetry parameter into n_values list */
    for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
        /* row of neighbor list */
        fNeighbors.RowAlias(i, neighbors);

        /* tags */
        int   tag_i = neighbors[0]; /* self is 1st spot */
        int  type_i = fType[tag_i];
        int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

        int valuep = 0;
        for (int is = 0; is < num_stresses; is++)
        {
            n_values(local_i,ndof+2+num_stresses+valuep++) = s_values(local_i,is);
        }

        /* recover J, the determinant of the deformation gradient, for atom i
		 * and divide stress values by it */
		double J = s_values(local_i,num_stresses);
		for (int is = 0; is < num_stresses; is++) 
			n_values(local_i,ndof+2+is) /= J;

        for (int n = 0; n < ndof; n++)
            n_values(local_i, ndof+2+num_stresses+num_stresses+n) = s_values(local_i,num_stresses+1+n);

        n_values(local_i, num_output-1) = s_values(local_i,num_s_vals-1);
    }

#endif /* NO_PARTICLE_STRESS_OUTPUT */

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

/* describe the parameters needed by the interface */
void EAMT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParticleT::DefineParameters(list);
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

/***********************************************************************
 * Protected
 ***********************************************************************/

/* generate labels for output data */
void EAMT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
  int ndof=NumDOF();
  if (ndof > 3) ExceptionT::GeneralFail("EAMT::GenerateOutputLabels");

  /* displacement labels */
  const char* disp[3] = {"D_X", "D_Y", "D_Z"};
  const char* SV[3] = {"SV_X", "SV_Y", "SV_Z"};
  int num_labels =
    ndof // displacements
    + 2;     // PE and KE

#ifndef NO_PARTICLE_STRESS_OUTPUT
	int num_stress=0;
	const char* stress[6];
	const char* strain[6];
	if (ndof==3){
	  num_stress=6;
	  stress[0]="s11";
	  stress[1]="s22";
	  stress[2]="s33";
	  stress[3]="s23";
	  stress[4]="s13";
	  stress[5]="s12";
	  }
	  else if (ndof==2) {
	   num_stress=3;
	  stress[0]="s11";
	  stress[1]="s22";
	  stress[2]="s12";
	  }
	  else if (ndof==1) {
	   num_stress=1;
	  stress[0] = "s11";
	  }
	if (ndof==3){
	  
	  strain[0]="e11";
	  strain[1]="e12";
	  strain[2]="e13";
	  strain[3]="e22";
	  strain[4]="e23";
	  strain[5]="e33";
	  }
	  else if (ndof==2) {
	   
	  strain[0]="e11";
	  strain[1]="e12";
	  strain[2]="e22";
	  }
	  else if (ndof==1) {

	  strain[0] = "e11";
	  }
	num_labels+=num_stress;
	num_labels++; //another label for the centrosymmetry
	num_labels+=num_stress; //another for the strain
	num_labels+=ndof; /*and another for the slip vector*/
#endif /* NO_PARTICLE_STRESS_OUTPUT */

  labels.Dimension(num_labels);
  int dex = 0;
  for (dex = 0; dex < NumDOF(); dex++)
    labels[dex] = disp[dex];
  labels[dex++] = "PE";
  labels[dex++] = "KE";

#ifndef NO_PARTICLE_STRESS_OUTPUT
	for (int ns =0 ; ns<num_stress; ns++)
	  labels[dex++]=stress[ns];
	for (int ns =0 ; ns<num_stress; ns++)
	  labels[dex++]=strain[ns];
	for (int i=0; i<ndof; i++)
	  labels[dex++]=SV[i];
	labels[dex++]= "CS";
#endif /* NO_PARTICLE_STRESS_OUTPUT */
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
  if (formM) 
    {
    /* collect mass per particle */
    dArrayT mass(fNumTypes);
    for (int i = 0; i < fNumTypes; i++)
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
  GenerateNeighborList(part_nodes, NearestNeighborDistance, NearestNeighbors, true, true);
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
  if (ElementSupport().PrintInput())
    {
      out << " Neighbor lists (self as leading neighbor):\n";
      out << setw(kIntWidth) << "row" << "  n..." << '\n';
      iArrayT tmp(fNeighbors.Length(), fNeighbors.Pointer());
      tmp++;
      fNeighbors.WriteNumbered(out);
      tmp--;
      out.flush();
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

/* construct the list of properties from the given input stream */
void EAMT::EchoProperties(ifstreamT& in, ofstreamT& out)
{
  /* read potentials : one potential corresponds to one type, 
                       cannot mix different potentials here*/
  int num_potentials = -1;
  in >> num_potentials;  

  fEAMProperties.Dimension(num_potentials); 
  fEAMProperties = NULL;
  for (int i = 0; i < fEAMProperties.Length(); i++)
    {
      int type_of_file;
      in >> type_of_file;

      StringT file;
      in >> file;
      file.ToNativePathName();
      
      StringT path;
      path.FilePath(in.filename());	
      file.Prepend(path);
      
      fEAMProperties[i] = new ParadynEAMT(file);
    }


  /* echo particle properties */
  out << "\n Particle properties:\n\n";
  out << " Number of properties. . . . . . . . . . . . . . = " 
      << fEAMProperties.Length() << '\n';
  for (int i = 0; i < fEAMProperties.Length(); i++)
    {
      out << " Property: " << i+1 << '\n';
      fEAMProperties[i]->Write(out);
    }
	
  /* copy into base class list */
  fParticleProperties.Dimension(fEAMProperties.Length());
  for (int i = 0; i < fEAMProperties.Length(); i++)
    fParticleProperties[i] = fEAMProperties[i];
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
