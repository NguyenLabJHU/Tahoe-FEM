/* $Id */
#include "EAMT.h"

#include "fstreamT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"
#include "dSPMatrixT.h"

/* EAM potentials */
#include "ParadynEAMT.h"

using namespace Tahoe;

/* parameters */
const int kMemoryHeadRoom = 15; /* percent */

/* constructor */
EAMT::EAMT(const ElementSupportT& support, const FieldT& field):
  ParticleT(support, field),
  fNeighbors(kMemoryHeadRoom),
  fEqnos(kMemoryHeadRoom),
  fForce_list_man(0, fForce_list)
{

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
  /* inherited */
  ParticleT::Initialize();

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

  /* number of nodes */
  const ArrayT<int>* parition_nodes = fCommManager.PartitionNodes();
  int non = (parition_nodes) ? parition_nodes->Length() : ElementSupport().NumNodes();

  /* map from partition node index */
  const InverseMapT* inverse_map = fCommManager.PartitionNodes_inv();

  /* output arrays length number of active nodes */
  dArray2DT n_values(non, num_output), e_values;
  n_values = 0.0;

  /* global coordinates */
  const dArray2DT& coords = ElementSupport().CurrentCoordinates();

  /* pair properties function pointers */
  int current_property_i = -1;
  int current_property_j = -1;
  EAMPropertyT::PairEnergyFunction  pair_energy = NULL;

  EAMPropertyT::EmbedEnergyFunction embed_energy_i= NULL;
  EAMPropertyT::EDEnergyFunction ed_energy_j   = NULL;

  /* the field */
  const FieldT& field = Field();
  const dArray2DT& displacement = field[0];
  const dArray2DT* velocities = NULL;
  if (field.Order() > 0) velocities = &(field[1]);

  /* collect displacements */
  dArrayT vec, values_i;
  for (int i = 0; i < non; i++) {
    int   tag_i = (parition_nodes) ? (*parition_nodes)[i] : i;
    int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

    /* values for particle i */
    n_values.RowAlias(local_i, values_i);

    /* copy in */
    vec.Set(ndof, values_i.Pointer());
    displacement.RowCopy(tag_i, vec);
  }

  /* collect mass per particle */
  dArrayT mass(fNumTypes);
  for (int i = 0; i < fNumTypes; i++)
    mass[i] = fEAMProperties[fPropertiesMap(i,i)]->Mass();
	
  /* run through neighbor list */
  iArrayT neighbors;
  dArrayT x_i, x_j, r_ij(ndof);
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* tags */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];		
      int local_i = (inverse_map) ? inverse_map->Map(tag_i) : tag_i;

      /* values for particle i */
      n_values.RowAlias(local_i, values_i);		

      /* kinetic energy */
      if (velocities)
	{
	  velocities->RowAlias(tag_i, vec);
	  values_i[ndof+1] = 0.5*mass[type_i]*dArrayT::Dot(vec, vec);
	}
      coords.RowAlias(tag_i, x_i);

      /* Get rho_i = \sum_{j=1,neighbor(i)} rho_j(r_ij)  */
      double rho_i = 0.0;
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* tags */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];		
			
	  /* set electron density (if not already set) */
	  int property_j = fPropertiesMap(type_j, type_i);  // rho_j <=> potential j, type j
	  if (property_j != current_property_j)
	    {
	      ed_energy_j = fEAMProperties[property_j]->getElecDensEnergy();
	      current_property_j = property_j;
	    }
	  
	  /* global coordinates */
	  coords.RowAlias(tag_j, x_j);
	  
	  /* connecting vector */
	  r_ij.DiffOf(x_j, x_i);
	  double r = r_ij.Magnitude();
	  
	  rho_i += ed_energy_j(r,NULL,NULL);
	}

      /* Get embedding energy */
      double E_i = 0.0;
      int property_i = fPropertiesMap(type_i, type_i);  // same type to get embedding energy : E_i(rho_i)
      embed_energy_i = fEAMProperties[property_i]->getEmbedEnergy();
      E_i = embed_energy_i(rho_i,NULL,NULL); 
      // do not save current_property because needs to be used with pair.

      values_i[ndof] = E_i;
      n_values(local_i, ndof) = E_i;

      /* Get pair energy */
      /* run though neighbors for one atom - first neighbor is self
       * to compute potential energy */
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* tags */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];		
			
	  /* set pair property (if not already set) */
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_energy  = fEAMProperties[property_i]->getPairEnergy();
	      current_property_i = property_i;
	    }
	  
	  /* global coordinates */
	  coords.RowAlias(tag_j, x_j);
	  
	  /* connecting vector */
	  r_ij.DiffOf(x_j, x_i);
	  double r = r_ij.Magnitude();
	  
	  /* split interaction energy */
	  double uby2 = 0.5*pair_energy(r, NULL, NULL);
	  values_i[ndof] += uby2;
	  
	  /* second node may not be on processor */
	  if (!proc_map || (*proc_map)[tag_j] == rank) 
	    {
	      int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;
	      
	      if (local_j < 0 || local_j >= n_values.MajorDim())
		cout << caller << ": out of range: " << local_j << '\n';
	      else
		n_values(local_j, ndof) += uby2;
	    }
	}
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

  /* pair properties function pointers */
  int current_property = -1;
  EAMPropertyT::PairForceFunction pair_force = NULL;
  EAMPropertyT::PairStiffnessFunction pair_stiffness = NULL;

  /* work space */
  dArrayT r_ij(NumDOF(), fRHS.Pointer());
  dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

  /* run through neighbor list */
  const iArray2DT& field_eqnos = Field().Equations();
  iArrayT row_eqnos, col_eqnos; 
  iArrayT neighbors;
  dArrayT x_i, x_j;
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* type */
      int  tag_i = neighbors[0]; /* self is 1st spot */
      int type_i = fType[tag_i];
		
      /* particle equations */
      field_eqnos.RowAlias(tag_i, row_eqnos);

      /* run though neighbors for one atom - first neighbor is self */
      coords.RowAlias(tag_i, x_i);
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int tag_j = neighbors[j];
			
	  /* particle is a target column */
	  int col_eq_index = col_to_col_eq_row_map.Map(tag_j);
	  if (col_eq_index != -1)
	    {
	      /* more particle info */
	      int type_j = fType[tag_j];

	      /* particle equations */
	      col_eq.RowAlias(col_eq_index, col_eqnos);

	      /* set pair property (if not already set) */
	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  pair_force = fEAMProperties[property]->getPairForce();
		  pair_stiffness = fEAMProperties[property]->getPairStiffness();
		  current_property = property;
		}

	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);

	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
	      r_ji.SetToScaled(-1.0, r_ij);

	      /* interaction functions */
	      double F = pair_force(r, NULL, NULL);
	      double K = pair_stiffness(r, NULL, NULL);
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
	}
    }
}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* generate labels for output data */
void EAMT::GenerateOutputLabels(ArrayT<StringT>& labels) const
{
  if (NumDOF() > 3) ExceptionT::GeneralFail("EAMT::GenerateOutputLabels");

  /* displacement labels */
  const char* disp[3] = {"D_X", "D_Y", "D_Z"};
	
  int num_labels =
    NumDOF() // displacements
    + 2;     // PE and KE

  labels.Dimension(num_labels);
  int dex = 0;
  for (dex = 0; dex < NumDOF(); dex++)
    labels[dex] = disp[dex];
  labels[dex++] = "PE";
  labels[dex++] = "KE";
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
    for (int type_i = 0; type_i < fNumTypes; type_i++)
      {
	int property = fPropertiesMap(type_i,type_i);
	mass[type_i] = fEAMProperties[property]->Mass();
      }
    mass *= constM;
	
    AssembleParticleMass(mass);
  }
	
  /* assemble diagonal stiffness */
  if (formK && sys_type == GlobalT::kDiagonal)
    {
      /* assembly information */
      const ElementSupportT& support = ElementSupport();
      int group = Group();
      int ndof = NumDOF();
	
      /* global coordinates */
      const dArray2DT& coords = support.CurrentCoordinates();

      /* EAM properties function pointers */
      int current_property = -1;
      EAMPropertyT::PairForceFunction pair_force = NULL;
      EAMPropertyT::PairStiffnessFunction pair_stiffness = NULL;

      EAMPropertyT::EmbedForceFunction embed_force = NULL;
      EAMPropertyT::EmbedStiffnessFunction embed_stiffness = NULL;

      EAMPropertyT::EDForceFunction ed_energy = NULL;
      EAMPropertyT::EDForceFunction ed_force = NULL;    
      EAMPropertyT::EDStiffnessFunction ed_stiffness = NULL;

      fForce = 0.0;
      iArrayT neighbors;
      /* Get the \sum_{j=1,neighbor(i)} rho_j (r_ij) */
      dArrayT sum_rho(fNeighbors.MajorDim());
      sum_rho = 0.0;

      /* Get the \sum_{j=1,neighbor(i)} rho_j (r_ij) */
      for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
	  /* row of neighbor list for i */
	  fNeighbors.RowAlias(i, neighbors);
	  /* type */
	  int   tag_i = neighbors[0]; /* self is 1st spot */
	  int  type_i = fType[tag_i];
	  double* f_i = fForce(tag_i);
	  double* x_i = coords(tag_i);

	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int   tag_j = neighbors[j];
	      int  type_j = fType[tag_j];
	      double* x_j = coords(tag_j);
	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  ed_energy = fEAMProperties[property]->getElecDensEnergy();
		  current_property = property;
		}
	      
	      /* connecting vector */
	      double r_ij_0 = x_j[0] - x_i[0];
	      double r_ij_1 = x_j[1] - x_i[1];
	      double r_ij   = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);
	      
	      sum_rho[i] += ed_energy(r_ij,NULL,NULL); // += : sum over j
	    }
	}


      fForce = 0.0;
      /* run through neighbor list */
      dArrayT x_i, x_j, r_ij(ndof);

      /* Get the \sum_{j=1,neighbor(i)} rho'_j(r_ij)  r_ij/r */
      dArray2DT sum_rhop_r(fNeighbors.MajorDim(),ndof);
      sum_rhop_r = 0.0;

      for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
	  /* row of neighbor list */
	  fNeighbors.RowAlias(i, neighbors);

	  /* type */
	  int  tag_i = neighbors[0]; /* self is 1st spot */
	  int type_i = fType[tag_i];
	  double* k_i = fForce(tag_i);
		
	  coords.RowAlias(tag_i, x_i);

	  dArrayT Delec_dens_j(neighbors.Length());
	  Delec_dens_j = 0.0;

	  dArrayT DDelec_dens_j(neighbors.Length());
	  DDelec_dens_j = 0.0;

	  /* Get electron density rho'(j) and hessian rho''(j) */
	  /* Get the \sum_{j=1,neighbor(i)} rho'_j(r_ij)  r_ij/r */
	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int   tag_j = neighbors[j];
	      int  type_j = fType[tag_j];
	      
	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  ed_force         = fEAMProperties[property]->getElecDensForce();
		  ed_stiffness     = fEAMProperties[property]->getElecDensStiffness();

		  current_property = property;
		}
	      
	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);

	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();

	      Delec_dens_j[j]  = ed_force(r,NULL,NULL);      // [rho_j(r_ij)]'
	      DDelec_dens_j[j] = ed_stiffness(r,NULL,NULL);  // [rho_j(r_ij)]''

	      for (int k = 0; k < ndof; k++)
		sum_rhop_r(i,k) += Delec_dens_j[j]*r_ij[k]/r; // += : sum over j
	    }


	  /* Get contribution to the stiffness from embedding and pair forces */
	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      double* k_j = fForce(tag_j);
			
	      /* set pair property (if not already set) */
	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  pair_force     = fEAMProperties[property]->getPairForce();
		  pair_stiffness = fEAMProperties[property]->getPairStiffness();

		  embed_force = fEAMProperties[property]->getEmbedForce();
		  embed_stiffness = fEAMProperties[property]->getEmbedStiffness();

		  ed_force    = fEAMProperties[property]->getElecDensForce();
		  ed_stiffness= fEAMProperties[property]->getElecDensStiffness();

		  current_property = property;
		}
		
	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);
		
	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
			
	      /* pair interaction functions */
	      double F = pair_force(r, NULL, NULL);
	      double K = pair_stiffness(r, NULL, NULL);
	      K = (K < 0.0) ? 0.0 : K;

	      double Fbyr = F/r;
	      for (int k = 0; k < ndof; k++)
		{
		  double r_k = r_ij[k]*r_ij[k]/r/r;
		  double K_k = constK*(K*r_k + Fbyr*(1.0 - r_k));
		  k_i[k] += K_k;
		  k_j[k] += K_k;
		}

	      /* embedding interaction functions */
	      double F_rho_i = embed_force(sum_rho[i],NULL,NULL); // [F_i(\sum_k rho(r_ik))]'
	      double F_rho_j = embed_force(sum_rho[j],NULL,NULL); // [F_j(\sum_k rho(r_kj))]'

	      double K_rho_i = embed_stiffness(sum_rho[i],NULL,NULL); // [F_i(\sum_k rho(r_ik))]''
	      double K_rho_j = embed_stiffness(sum_rho[j],NULL,NULL); // [F_j(\sum_k rho(r_kj))]''
		
	      double Delec_dens_i = ed_force(r,NULL,NULL);
	      F = Delec_dens_i * F_rho_i + Delec_dens_j[j] * F_rho_j;
	      Fbyr = F/r;

	      double DDelec_dens_i = ed_stiffness(r,NULL,NULL);
	      K  = DDelec_dens_i * F_rho_j + DDelec_dens_j[j] * F_rho_i +
		   Delec_dens_i * Delec_dens_i * K_rho_j;
	      double K2 = DDelec_dens_j[j] * K_rho_i ;

	      dArrayT K2_k(ndof);
	      for (int k = 0; k < ndof; k++)
		  K2_k[k] = constK*(K2 * sum_rhop_r(i,k)*r_ij[k]/r);
	      	      
	      for (int k = 0; k < ndof; k++)
		{
		  double r_k = r_ij[k]*r_ij[k]/r/r;
		  double K_k = constK*(K*r_k + Fbyr*(1.0 - r_k));
		  k_i[k] += K_k + K2_k[k];   // += : sum over j
		  k_j[k] += K_k + K2_k[k];
		}
	    }
	}

      /* assemble */
      support.AssembleLHS(group, fForce, Field().Equations());
    }
  else if (formK)
    {
      /* assembly information */
      const ElementSupportT& support = ElementSupport();
      int group = Group();
      int ndof = NumDOF();
      fLHS.Dimension(2*ndof);
		
      /* global coordinates */
      const dArray2DT& coords = support.CurrentCoordinates();

      /* EAM properties function pointers */
      int current_property = -1;
      EAMPropertyT::PairForceFunction pair_force = NULL;
      EAMPropertyT::PairStiffnessFunction pair_stiffness = NULL;

      EAMPropertyT::EmbedForceFunction embed_force = NULL;
      EAMPropertyT::EmbedStiffnessFunction embed_stiffness = NULL;

      EAMPropertyT::EDForceFunction ed_energy = NULL;
      EAMPropertyT::EDForceFunction ed_force = NULL;    
      EAMPropertyT::EDStiffnessFunction ed_stiffness = NULL;

      /* work space */
      dArrayT r_ij(NumDOF(), fRHS.Pointer());
      dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

      /* run through neighbor list */
      const iArray2DT& field_eqnos = Field().Equations();
      iArray2DT pair_eqnos(2, ndof); 
      iArrayT pair(2);
      iArrayT neighbors;
      dArrayT x_i, x_j;

      /* Get the \sum_{j=1,neighbor(i)} rho_j (r_ij) */
      dArrayT sum_rho(fNeighbors.MajorDim());
      sum_rho = 0.0;
      for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
	  /* row of neighbor list */
	  fNeighbors.RowAlias(i, neighbors);

	  /* type */
	  int  tag_i = neighbors[0]; /* self is 1st spot */
	  int type_i = fType[tag_i];
	  pair[0] = tag_i;
		
	  coords.RowAlias(tag_i, x_i);
	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      pair[1] = tag_j;
			
	      /* set pair property (if not already set) */
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

	      sum_rho[i] += ed_energy(r,NULL,NULL);	      
	    }
	}

      /* Get the \sum_{j=1,neighbor(i)} rho'_j(r_ij)  r_ij/r */
      dArray2DT sum_rhop_r(fNeighbors.MajorDim(),ndof);
      sum_rhop_r = 0.0;

      /* Get the \sum_{j=1,neighbor(i)} rho'_i(r_ij)  F''_j  r_ij/r */
      dArray2DT sum_Fpp_rhop_r(fNeighbors.MajorDim(),ndof);
      sum_Fpp_rhop_r = 0.0;

      for (int i = 0; i < fNeighbors.MajorDim(); i++)
	{
	  /* row of neighbor list */
	  fNeighbors.RowAlias(i, neighbors);

	  /* type */
	  int  tag_i = neighbors[0]; /* self is 1st spot */
	  int type_i = fType[tag_i];
	  pair[0] = tag_i;

	  coords.RowAlias(tag_i, x_i);
	  /* Get electron density derivative rho'(j) and hessian rho''(j) */
	  dArrayT Delec_dens_j(neighbors.Length());
	  Delec_dens_j = 0.0;

	  /* Get electron density rho'(j) and hessian rho''(j) */
	  /* Get the \sum_{j=1,neighbor(i)} rho'_j(r_ij)  r_ij/r */
	  /* Get the \sum_{j=1,neighbor(i)} rho'_i(r_ij)  F''_j  r_ij/r */
	  dArrayT DDelec_dens_j(neighbors.Length());
	  DDelec_dens_j = 0.0;

	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      pair[1] = tag_j;
			
	      /* set pair property (if not already set) */
	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  pair_force     = fEAMProperties[property]->getPairForce();
		  pair_stiffness = fEAMProperties[property]->getPairStiffness();

		  embed_force    = fEAMProperties[property]->getEmbedForce();
		  embed_stiffness= fEAMProperties[property]->getEmbedStiffness();

		  ed_force       = fEAMProperties[property]->getElecDensForce();
		  ed_stiffness   = fEAMProperties[property]->getElecDensStiffness();

		  current_property = property;
		}
		
	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);
		
	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
	      r_ji.SetToScaled(-1.0, r_ij);
			
	      Delec_dens_j[j]  = ed_force(r,NULL,NULL);
	      DDelec_dens_j[j] = ed_stiffness(r,NULL,NULL);

	      for (int k = 0; k < ndof; k++)
		{
		  sum_rhop_r(i,k) += Delec_dens_j[j]*r_ij[k]/r; // += : sum over j
		  //sum_Fpp_rhop_r(i,k) += 
		}
	    }


	  /* Get contribution to the stiffness from embedding and pair forces */
	  for (int j = 1; j < neighbors.Length(); j++)
	    {
	      /* global tag */
	      int  tag_j = neighbors[j];
	      int type_j = fType[tag_j];
	      pair[1] = tag_j;
			
	      /* set pair property (if not already set) */
	      int property = fPropertiesMap(type_i, type_j);
	      if (property != current_property)
		{
		  pair_force     = fEAMProperties[property]->getPairForce();
		  pair_stiffness = fEAMProperties[property]->getPairStiffness();

		  embed_force    = fEAMProperties[property]->getEmbedForce();
		  embed_stiffness= fEAMProperties[property]->getEmbedStiffness();

		  ed_force       = fEAMProperties[property]->getElecDensForce();
		  ed_stiffness   = fEAMProperties[property]->getElecDensStiffness();

		  current_property = property;
		}
		
	      /* global coordinates */
	      coords.RowAlias(tag_j, x_j);
		
	      /* connecting vector */
	      r_ij.DiffOf(x_j, x_i);
	      double r = r_ij.Magnitude();
	      r_ji.SetToScaled(-1.0, r_ij);
			
	      /* pair interaction functions */
	      double F = constK*pair_force(r, NULL, NULL);
	      double K = constK*pair_stiffness(r, NULL, NULL);
	      double Fbyr = F/r;

	      /* 1st term */
	      fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);
		
	      /* 2nd term */
	      fLHS.AddScaled(Fbyr, fOneOne);
				
	      /* embedding interaction functions */
	      double F_rho_i = embed_force(sum_rho[j],NULL,NULL);
	      double F_rho_j = embed_force(sum_rho[i],NULL,NULL);

	      double K_rho_i = embed_stiffness(sum_rho[i],NULL,NULL);
	      double K_rho_j = embed_stiffness(sum_rho[j],NULL,NULL);

	      double Delec_dens_i = ed_force(r,NULL,NULL);
	      F = Delec_dens_i * F_rho_i + Delec_dens_j[j] * F_rho_j;
	      double DDelec_dens_i = ed_stiffness(r,NULL,NULL);
	      K = DDelec_dens_j[j] * F_rho_i + DDelec_dens_i * F_rho_j;
	      Fbyr = F/r;

	      /* 1st term */
	      fLHS.Outer(fRHS, fRHS, (K - Fbyr)/r/r);
		
	      /* 2nd term */
	      fLHS.AddScaled(Fbyr, fOneOne);

	      /* last term */
	      for (int i = 0; i < ndof; i++)
		for (int j = i+1; j < ndof; j++)
		  {

		    for (int k = 0; k < ndof; k++)
		      fLHS(i,j) += 0.0;

		  }


				
	      /* assemble */
	      pair_eqnos.RowCollect(pair, field_eqnos);
	      support.AssembleLHS(group, fLHS, pair_eqnos);


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

  /* EAM properties function pointers */
  int current_property_i = -1;
  int current_property_j = -1;
  EAMPropertyT::EDEnergyFunction ed_energy_i = NULL,ed_energy_j = NULL;
  EAMPropertyT::EDForceFunction  ed_force_i = NULL,ed_force_j = NULL;
	
  fForce = 0.0;
  iArrayT neighbors;

  /* Get rho_i = \sum_{k=1,neighbor(i)} rho_k (r_ik)  */
  dArrayT rho_i(fNeighbors.MajorDim());
  rho_i = 0.0;
  /* Get rho_j = \sum_{k=1,neighbor(j)} rho_k (r_jk)  */
  dArrayT rho_j(fNeighbors.MajorDim());
  rho_j = 0.0;
  /* run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);
      
      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      double* f_i = fForce(tag_i);
      double* x_i = coords(tag_i);
		
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  double* x_j = coords(tag_j);

	  int property_i = fPropertiesMap(type_i, type_j);  // type_i <=> rho_i
	  if (property_i != current_property_i)
	    {
	      ed_energy_i  = fEAMProperties[property_i]->getElecDensForce();
	      current_property_i = property_i;
	    }

	  int property_j = fPropertiesMap(type_j, type_i);  // type_j <=> rho_j
	  if (property_j != current_property_j)
	    {
	      ed_energy_j  = fEAMProperties[property_j]->getElecDensForce();
	      current_property_j = property_j;
	    }
		
	  /* connecting vector */
	  double r_ij_0 = x_j[0] - x_i[0];
	  double r_ij_1 = x_j[1] - x_i[1];
	  double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);
		
	  rho_i[i] += ed_energy_i(r,NULL,NULL); // sum over j
	  rho_j[j] += ed_energy_j(r,NULL,NULL); // sum over i
	}
    }
  current_property_i = -1;
  current_property_j = -1;

  EAMPropertyT::PairForceFunction pair_force = NULL;
  EAMPropertyT::EmbedForceFunction embed_force_i = NULL,embed_force_j=NULL;

  fForce = 0.0;
  /* Compute Force : sum up contributions  */
  /* run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);
      
      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      double* f_i = fForce(tag_i);
      double* x_i = coords(tag_i);

      /* Compute Force  */
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  double* f_j = fForce(tag_j);
	  double* x_j = coords(tag_j);

	  /* set EAM property (if not already set) */
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_force    = fEAMProperties[property_i]->getPairForce();
	      embed_force_i = fEAMProperties[property_i]->getEmbedForce();
	      ed_force_i    = fEAMProperties[property_i]->getElecDensForce();

	      current_property_i = property_i;
	    }

	  int property_j = fPropertiesMap(type_j, type_i);
	  if (property_j != current_property_j)
	    {
	      embed_force_j = fEAMProperties[property_j]->getEmbedForce();
	      ed_force_j   = fEAMProperties[property_j]->getElecDensForce();

	      current_property_j = property_j;
	    }

	  /* connecting vector */
	  double r_ij_0 = x_j[0] - x_i[0];
	  double r_ij_1 = x_j[1] - x_i[1];
	  double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1);
			
	  /* Component of force coming from Pair potential */
	  double F = pair_force(r, NULL, NULL);
	  double Fbyr = formKd*F/r;

	  r_ij_0 *= Fbyr;
	  f_i[0] += r_ij_0;
	  f_j[0] +=-r_ij_0;

	  r_ij_1 *= Fbyr;
	  f_i[1] += r_ij_1;
	  f_j[1] +=-r_ij_1;

	  /* Component of force coming from Embedding energy */
	  double Ep_i   = embed_force_i(rho_i[i],NULL,NULL);
	  double rhop_i = ed_force_i(r,NULL,NULL);
	  double Ep_j   = embed_force_j(rho_j[j],NULL,NULL);
	  double rhop_j = ed_force_j(r,NULL,NULL);

	  F = rhop_i * Ep_j + rhop_i * Ep_i;
	  Fbyr = formKd*F/r;

	  r_ij_0 *= Fbyr;
	  f_i[0] += r_ij_0;
	  f_j[0] +=-r_ij_0;

	  r_ij_1 *= Fbyr;
	  f_i[1] += r_ij_1;
	  f_j[1] +=-r_ij_1;
	}
    }

    /* assemble */
    support.AssembleRHS(group, fForce, Field().Equations());
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

  /* pair properties function pointers */
  int current_property_i = -1;
  int current_property_j = -1;
  EAMPropertyT::EDEnergyFunction ed_energy_i = NULL,ed_energy_j = NULL;
  EAMPropertyT::EDForceFunction  ed_force_i = NULL,ed_force_j = NULL;

  fForce = 0.0;
  iArrayT neighbors;

  /* Get rho_i = \sum_{k=1,neighbor(i)} rho_k (r_ik)  */
  dArrayT rho_i(fNeighbors.MajorDim());
  rho_i = 0.0;
  /* Get rho_j = \sum_{k=1,neighbor(j)} rho_k (r_jk)  */
  dArrayT rho_j(fNeighbors.MajorDim());
  rho_j = 0.0;
  /* run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list for i */
      fNeighbors.RowAlias(i, neighbors);

      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      double* f_i = fForce(tag_i);
      double* x_i = coords(tag_i);

      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  double* x_j = coords(tag_j);

	  int property_i = fPropertiesMap(type_i, type_j);  // type_i <=> rho_i
	  if (property_i != current_property_i)
	    {
	      ed_energy_i  = fEAMProperties[property_i]->getElecDensForce();
	      current_property_i = property_i;
	    }

	  int property_j = fPropertiesMap(type_j, type_i);  // type_j <=> rho_j
	  if (property_j != current_property_j)
	    {
	      ed_energy_j  = fEAMProperties[property_j]->getElecDensForce();
	      current_property_j = property_j;
	    }

	  /* connecting vector */
	  double r_ij_0 = x_j[0] - x_i[0];
	  double r_ij_1 = x_j[1] - x_i[1];
	  double r_ij_2 = x_j[2] - x_i[2];
	  double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1 + r_ij_2*r_ij_2);
	  
	  rho_i[i] += ed_energy_i(r,NULL,NULL); // sum over j
	  rho_j[j] += ed_energy_j(r,NULL,NULL); // sum over i
	}
    }
  current_property_i = -1;
  current_property_j = -1;

  EAMPropertyT::PairForceFunction pair_force = NULL;
  EAMPropertyT::EmbedForceFunction embed_force_i = NULL,embed_force_j=NULL;

  fForce = 0.0;
  /* Compute Force : sum up contributions  */
  /* run through neighbor list */
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      /* row of neighbor list */
      fNeighbors.RowAlias(i, neighbors);

      /* type */
      int   tag_i = neighbors[0]; /* self is 1st spot */
      int  type_i = fType[tag_i];
      double* f_i = fForce(tag_i);
      double* x_i = coords(tag_i);
		
      /* Compute Force  */
      for (int j = 1; j < neighbors.Length(); j++)
	{
	  /* global tag */
	  int   tag_j = neighbors[j];
	  int  type_j = fType[tag_j];
	  double* f_j = fForce(tag_j);
	  double* x_j = coords(tag_j);

	  /* set EAM property (if not already set) */
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_force    = fEAMProperties[property_i]->getPairForce();
	      embed_force_i = fEAMProperties[property_i]->getEmbedForce();
	      ed_force_i    = fEAMProperties[property_i]->getElecDensForce();

	      current_property_i = property_i;
	    }

	  int property_j = fPropertiesMap(type_j, type_i);
	  if (property_j != current_property_j)
	    {
	      embed_force_j = fEAMProperties[property_j]->getEmbedForce();
	      ed_force_j   = fEAMProperties[property_j]->getElecDensForce();

	      current_property_j = property_j;
	    }
		
	  /* connecting vector */
	  double r_ij_0 = x_j[0] - x_i[0];
	  double r_ij_1 = x_j[1] - x_i[1];
	  double r_ij_2 = x_j[2] - x_i[2];
	  double r      = sqrt(r_ij_0*r_ij_0 + r_ij_1*r_ij_1 + r_ij_2*r_ij_2);
			
	  /* Component of force coming from Pair potential */
	  double F = pair_force(r, NULL, NULL);
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

	  /* Component of force coming from Embedding energy */
	  double Ep_i   = embed_force_i(rho_i[i],NULL,NULL);
	  double rhop_i = ed_force_i(r,NULL,NULL);
	  double Ep_j   = embed_force_j(rho_j[j],NULL,NULL);
	  double rhop_j = ed_force_j(r,NULL,NULL);

	  F = rhop_i * Ep_j + rhop_i * Ep_i;
	  Fbyr = formKd*F/r;

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
    }

  /* assemble */
  support.AssembleRHS(group, fForce, Field().Equations());
}

/* set neighborlists */
void EAMT::SetConfiguration(void)
{
  /* inherited */
  ParticleT::SetConfiguration();

  /* reset neighbor lists */
  const ArrayT<int>* part_nodes = fCommManager.PartitionNodes();
  if (fActiveParticles) 
    part_nodes = fActiveParticles;
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
}

/* construct the list of properties from the given input stream */
void EAMT::EchoProperties(ifstreamT& in, ofstreamT& out)
{
  /* read potentials : each potential corresponds here to a type ex: Ni or Cu*/
  int num_potentials = -1;
  in >> num_potentials;  

  fEAMProperties.Dimension(num_potentials); 
  fEAMProperties = NULL;
  for (int i = 0; i < fEAMProperties.Length(); i++)
    {
      ParticlePropertyT::TypeT property;
      in >> property;
      switch (property)
	{
	case ParticlePropertyT::kParadynEAM:
	  {
	    StringT file;
	    in >> file;
	    file.ToNativePathName();

	    StringT path;
	    path.FilePath(in.filename());	
	    file.Prepend(path);
	  
	    fEAMProperties[i] = new ParadynEAMT(file);
	    break;
	  }
	default:
	  ExceptionT::BadInputValue("EAMT::ReadProperties", 
				    "unrecognized property type: %d", 
                                     property);
	}
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
