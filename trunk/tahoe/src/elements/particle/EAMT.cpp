/* $Id: EAMT.cpp,v 1.37 2003-07-15 23:31:49 saubry Exp $ */
#include "EAMT.h"

#include "fstreamT.h"
#include "eIntegratorT.h"
#include "InverseMapT.h"
#include "CommManagerT.h"
#include "dSPMatrixT.h"

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
  fEqnos(kMemoryHeadRoom),
  fForce_list_man(0, fForce_list),
  fElectronDensity_man(kMemoryHeadRoom, fElectronDensity, 1),
  fEmbeddingEnergy_man(kMemoryHeadRoom, fEmbeddingEnergy, 1),
  fEmbeddingForce_man(kMemoryHeadRoom, fEmbeddingForce, 1),
  fEmbeddingStiff_man(kMemoryHeadRoom, fEmbeddingStiff, 1),
  frhop_r_man(kMemoryHeadRoom, frhop_r,NumDOF())
{}

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

  /* inherited */
  ParticleT::Initialize();

  /* dimension */
  int ndof = NumDOF();
  fLHS.Dimension(2*ndof);
  fRHS.Dimension(2*ndof);

  /* new array */
  fRHS2.Dimension(2*ndof);

  /* constant matrix needed to calculate stiffness */
  fOneOne.Dimension(fLHS);
  dMatrixT one(ndof);
  one.Identity();
  fOneOne.SetBlock(0, 0, one);
  fOneOne.SetBlock(ndof, ndof, one);
  one *= -1;
  fOneOne.SetBlock(0, ndof, one);
  fOneOne.SetBlock(ndof, 0, one);

  /* set up communication of electron density information */
  fElectronDensityMessageID = fCommManager.Init_AllGather(MessageT::Double, 1);
  fElectronDensity_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of embedding energy information */
  fEmbeddingEnergyMessageID = fCommManager.Init_AllGather(MessageT::Double, 1);
  fEmbeddingEnergy_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of embedding force information */
  fEmbeddingForceMessageID = fCommManager.Init_AllGather(MessageT::Double, 1);
  fEmbeddingForce_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of embedding stiffness information */
  fEmbeddingStiffMessageID = fCommManager.Init_AllGather(MessageT::Double, 1);
  fEmbeddingStiff_man.SetMajorDimension(ElementSupport().NumNodes(), false);

  /* set up communication of rhop * r information */
  frhop_rMessageID = fCommManager.Init_AllGather(MessageT::Double, NumDOF());
  frhop_r_man.SetMajorDimension(ElementSupport().NumNodes(), false);
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

  /* the field */
  const FieldT& field = Field();
  const dArray2DT& displacement = field[0];
  const dArray2DT* velocities = NULL;
  if (field.Order() > 0) velocities = &(field[1]);

  /* collect displacements */
  dArrayT vec, values_i;
  for (int i = 0; i < non; i++) 
    {
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

  if(iEmb == 1)
    {
      /* get electron density */
      if(ndof == 2) GetRho2D(coords,fElectronDensity);
      if(ndof == 3) GetRho3D(coords,fElectronDensity);
      /* exchange electron density information */
      fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);
      
      /* get embedding energy */
      GetEmbEnergy(coords,fElectronDensity,fEmbeddingEnergy);
      /* exchange embedding energy information */
      fCommManager.AllGather(fEmbeddingEnergyMessageID, fEmbeddingEnergy);
    }

  /* EAM properties function pointers */
  EAMPropertyT::PairEnergyFunction  pair_energy_i = NULL;
  EAMPropertyT::PairEnergyFunction  pair_energy_j = NULL;

  iArrayT neighbors;
  dArrayT x_i, x_j, r_ij(ndof);

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
			
	  int property_i = fPropertiesMap(type_i, type_j);
	  if (property_i != current_property_i)
	    {
	      pair_energy_i  = fEAMProperties[property_i]->getPairEnergy();
	      current_property_i = property_i;
	    }
	  
	  int property_j = fPropertiesMap(type_j, type_i);
	  if (property_j != current_property_j)
	    {
	      pair_energy_j  = fEAMProperties[property_j]->getPairEnergy();
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
	  
	  /* second node may not be on processor */
	  if (!proc_map || (*proc_map)[tag_j] == rank) 
	    {
	      int local_j = (inverse_map) ? inverse_map->Map(tag_j) : tag_j;
	      
	      if (local_j < 0 || local_j >= n_values.MajorDim())
		cout << caller << ": out of range: " << local_j << '\n';
	      else
		n_values(local_j, ndof) += phiby2;
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

  /* work space */
  dArrayT r_ij(NumDOF(), fRHS.Pointer());
  dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

  /* run through neighbor list */
  const iArray2DT& field_eqnos = Field().Equations();
  iArrayT row_eqnos, col_eqnos; 
  iArrayT neighbors;
  dArrayT x_i, x_j;

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
		  rp_i[k] += rhop_i * r_ij[k]; 
		  rp_j[k] += rhop_j * r_ij[k]; 
		}
	    }
	}
    }
 
  /* exchange electron density information */
  fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);

  /* exchange rhop * r information */
  fCommManager.AllGather(frhop_rMessageID, frhop_r);

  /* get embedding force */
  GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
  /* exchange embedding force information */
  fCommManager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

  /* get embedding stiffness */
  GetEmbStiff(coords,fElectronDensity,fEmbeddingStiff);
  /* exchange embedding stiffness information */
  fCommManager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);

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
		  
		  double E =  0.5* z_i * z_j/r;
		  double F =  0.5*(z_i * zp_j + zp_i * z_j)/r - E/r;
		  double K =  0.5*(zpp_i*z_j + 2*zp_i * zp_j + z_i * zpp_j)/r - 2*F/r;
		  
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
		  r_ji = r_ij;

		  double Ep_i   = fEmbeddingForce(tag_i,0); 
		  double Epp_i  = fEmbeddingStiff(tag_i,0); 
		  
		  double rhop_i  = ed_force_i(r,NULL,NULL);
		  double rhopp_i = ed_stiffness_i(r,NULL,NULL);
		  
		  double Ep_j   = fEmbeddingForce(tag_j,0); 
		  double Epp_j  = fEmbeddingStiff(tag_j,0); 
		  
		  double rhop_j  = ed_force_j(r,NULL,NULL);
		  double rhopp_j = ed_stiffness_j(r,NULL,NULL);
		  
		  double F_i = rhop_j * Ep_i;
		  double F_i_byr = F_i/r;

		  double F_j = rhop_i * Ep_j;
		  double F_j_byr = F_j/r;

		  double K_i = Ep_i * rhopp_j;
		  double K_j = Ep_j * rhopp_i;

		  double K2_i = Epp_i * rhop_j;
		  double K2_j = Epp_j * rhop_i;

		  /* 2nd term */
		  dArrayT El(2*ndof);
		  for (int k = 0; k < ndof; k++)
		    {
		      El[k]      =  (K_i - F_i_byr)*fRHS[k]/r + K2_i*rp_i[k];
		      El[k+ndof] =  (K_j - F_j_byr)*fRHS[k]/r + K2_j*rp_j[k];
		    }		  

		  fLHS.Outer(fRHS, El, 1./r);		  

		  /* 1st term */
		  for (int k = 0; k < ndof; k++)
		    {
		      fLHS(k,k)           +=  F_i_byr;
		      fLHS(k+ndof,k+ndof) +=  F_j_byr;
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
    for (int i = 0; i < fNumTypes; i++)
      mass[i] = fEAMProperties[fPropertiesMap(i,i)]->Mass();
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

      iArrayT neighbors;
      dArrayT x_i, x_j, r_ij(ndof);

      /* get electron density */
      if(ndof == 2) GetRho2D(coords,fElectronDensity);
      if(ndof == 3) GetRho3D(coords,fElectronDensity);
      /* exchange electron density information */
      fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);

      /* get embedding force */
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
      /* exchange embedding force information */
      fCommManager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

      /* get embedding stiffness */
      GetEmbStiff(coords,fElectronDensity,fEmbeddingStiff);
      /* exchange embedding stiffness information */
      fCommManager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);
   
      /* get rhop * r */
      frhop_r = 0.0;
      GetRhop_r(coords,frhop_r);
      /* exchange rhop * r information */
      fCommManager.AllGather(frhop_rMessageID, frhop_r);

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
		  
		  double E = 0.5*z_i * z_j/r;
		  double F = 0.5*(z_i * zp_j + zp_i * z_j)/r - E/r;
		  double K = 0.5*(zpp_i*z_j + 2*zp_i * zp_j + z_i * zpp_j)/r - 2*F/r;

		  double Fbyr = F/r;
		  
		  K = (K < 0.0) ? 0.0 : K;

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
		  
		  double F_i = rhop_j * Ep_i;
		  double F_i_byr = F_i/r;

		  double F_j = rhop_i * Ep_j;
		  double F_j_byr = F_j/r;

		  double K_i = Ep_i * rhopp_j;
		  double K_j = Ep_j * rhopp_i;

		  double K2_i = Epp_i * rhop_j;
		  double K2_j = Epp_j * rhop_i;

		  for (int k = 0; k < ndof; k++)
		    {
		      double r_k  = r_ij[k]*r_ij[k]/(r*r);
		      double r2_k = r_ij[k]/r;

		      double Fi = F_i_byr*(1.0 - r_k);
		      double Fj = F_j_byr*(1.0 - r_k);

		      double K_ik = K_i*r_k + K2_i*rp_i[k]*r2_k;
		      double K_jk = K_j*r_k + K2_j*rp_j[k]*r2_k;
		      
		      k_i[k] += constK*(Fi + K_ik);
		      k_j[k] += constK*(Fj + K_jk);
		    }
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

      /* work space */
      dArrayT r_ij(NumDOF(), fRHS.Pointer());
      dArrayT r_ji(NumDOF(), fRHS.Pointer() + NumDOF());

      const iArray2DT& field_eqnos = Field().Equations();
      iArray2DT pair_eqnos(2, ndof); 
      iArrayT pair(2);
      iArrayT neighbors;
      dArrayT x_i, x_j;

      /* component 11 */
      dArrayT D2fDx2_0(ElementSupport().NumNodes());
      D2fDx2_0 = 0.0;
      /* component 22 */
      dArrayT D2fDx2_1(ElementSupport().NumNodes());
      D2fDx2_1 = 0.0;
      /* component 33 */
      dArrayT D2fDx2_2(ElementSupport().NumNodes());
      D2fDx2_2 = 0.0;
      /* component 12 */
      dArrayT D2fDx2_3(ElementSupport().NumNodes());
      D2fDx2_3 = 0.0;
      /* component 13 */
      dArrayT D2fDx2_4(ElementSupport().NumNodes());
      D2fDx2_4 = 0.0;
      /* component 23 */
      dArrayT D2fDx2_5(ElementSupport().NumNodes());
      D2fDx2_5 = 0.0;

      /* EAM properties function pointers */
      if(ndof == 2) GetRho2D(coords,fElectronDensity);
      if(ndof == 3) GetRho3D(coords,fElectronDensity);
      /* exchange electron density information */
      fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);

      /* get embedding force */
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
      /* exchange embedding force information */
      fCommManager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

      /* get embedding stiffness */
      GetEmbStiff(coords,fElectronDensity,fEmbeddingStiff);
      /* exchange embedding stiffness information */
      fCommManager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);

      /* get rhop * r */
      frhop_r = 0.0;
      GetRhop_r(coords,frhop_r);
      /* exchange rhop * r information */
      fCommManager.AllGather(frhop_rMessageID, frhop_r);

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
		
		double E =  0.5* z_i * z_j/r;
		double F =  0.5*(z_i * zp_j + zp_i * z_j)/r - E/r;
		double K =  0.5*(zpp_i*z_j + 2*zp_i * zp_j + z_i * zpp_j)/r - 2*F/r;
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
		  r_ji = r_ij;

		  double Ep_i   = fEmbeddingForce(tag_i,0); 
		  double Epp_i  = fEmbeddingStiff(tag_i,0); 
		  
		  double rhop_i  = ed_force_i(r,NULL,NULL);
		  double rhopp_i = ed_stiffness_i(r,NULL,NULL);
		  
		  double Ep_j   = fEmbeddingForce(tag_j,0); 
		  double Epp_j  = fEmbeddingStiff(tag_j,0); 
		  
		  double rhop_j  = ed_force_j(r,NULL,NULL);
		  double rhopp_j = ed_stiffness_j(r,NULL,NULL);
		  
		  double F_i = rhop_j * Ep_i;
		  double F_i_byr = F_i/r;

		  double F_j = rhop_i * Ep_j;
		  double F_j_byr = F_j/r;

		  double K_i = Ep_i * rhopp_j;
		  double K_j = Ep_j * rhopp_i;

		  double K2_i = Epp_i * rhop_j;
		  double K2_j = Epp_j * rhop_i;

		  /* 2nd term */
		  dArrayT El(2*ndof);
		  for (int k = 0; k < ndof; k++)
		    {
		      El[k]      =  (K_i - F_i_byr)*fRHS[k]/r + K2_i*rp_i[k];
		      El[k+ndof] =  (K_j - F_j_byr)*fRHS[k]/r + K2_j*rp_j[k];
		    }		  

		  fLHS.Outer(fRHS, El, 1./r);		  

		  /* 1st term */
		  for (int k = 0; k < ndof; k++)
		    {
		      fLHS(k,k)           +=  F_i_byr;
		      fLHS(k+ndof,k+ndof) +=  F_j_byr;
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

  if(iEmb == 1)
    {
      /* get electron density */
      GetRho2D(coords,fElectronDensity);
      /* exchange electron density information */
      fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);
      
      /* get embedding force */
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);

      /* exchange embedding energy information */
      fCommManager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);
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
      double* x_i = coords(tag_i);

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
	      
	      double E = 0.5*z_i*z_j/r;
	      double F = 0.5*(z_i*zp_j + zp_i*z_j)/r - E/r;
	      
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
	      double rhop_i = ed_force_i(r,NULL,NULL);
	      
	      double Ep_j   = fEmbeddingForce(tag_j,0);
	      double rhop_j = ed_force_j(r,NULL,NULL);
	      
	      double F1 =  rhop_j * Ep_i;
	      double F1byr = formKd*F1/r;

	      double F2 =  rhop_i * Ep_j ;
	      double F2byr = formKd*F2/r;

	      f_i[0] += r_ij_0 * F1byr;
	      f_j[0] += r_ij_0 * F2byr;
	      
	      f_i[1] += r_ij_1 * F1byr;
	      f_j[1] += r_ij_1 * F2byr;
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
  
  if(iEmb == 1)
    {
      /* get electron density */
      fElectronDensity = 0.0;
      GetRho3D(coords,fElectronDensity);
      /* exchange electron density information */
      fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);
      
      /* get embedding force */
      fEmbeddingForce = 0.0;
      GetEmbForce(coords,fElectronDensity,fEmbeddingForce);
      /* exchange embedding energy information */
      fCommManager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);
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
	      
	      double E = 0.5*z_i*z_j/r;
	      double F = 0.5*(z_i*zp_j + zp_i*z_j)/r - E/r;
	      
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
	      double rhop_i = ed_force_i(r,NULL,NULL);
	      
	      double Ep_j   = fEmbeddingForce(tag_j,0);
	      double rhop_j = ed_force_j(r,NULL,NULL);
	      
	      double F1 =  rhop_j * Ep_i;
	      double F1byr = formKd*F1/r;

	      double F2 =  rhop_i * Ep_j ;
	      double F2byr = formKd*F2/r;

	      f_i[0] +=  r_ij_0 * F1byr;
	      f_j[0] += -r_ij_0 * F2byr;
	      
	      f_i[1] +=  r_ij_1 * F1byr;
	      f_j[1] += -r_ij_1 * F2byr;

	      f_i[2] +=  r_ij_2 * F1byr;
	      f_j[2] += -r_ij_2 * F2byr;	  
	  
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

  int nnd = ElementSupport().NumNodes();
  // ELECTRON DENSITY //
  /* reset the electron density array */
  fElectronDensity_man.SetMajorDimension(nnd, true);
  fElectronDensity.Resize(nnd);
  /* exchange type information */
  fCommManager.AllGather(fElectronDensityMessageID, fElectronDensity);

  // EMBEDDING ENERGY //
  /* reset the embedding energy array */
  fEmbeddingEnergy_man.SetMajorDimension(nnd, true);
  fEmbeddingEnergy.Resize(nnd);
  /* exchange type information */
  fCommManager.AllGather(fEmbeddingEnergyMessageID, fEmbeddingEnergy);

  /* reset the embedding force array */
  fEmbeddingForce_man.SetMajorDimension(nnd, true);
  fEmbeddingForce.Resize(nnd);
  /* exchange type information */
  fCommManager.AllGather(fEmbeddingForceMessageID, fEmbeddingForce);

  /* reset the embedding stiffness array */
  fEmbeddingStiff_man.SetMajorDimension(nnd, true);
  fEmbeddingStiff.Resize(nnd);
  /* exchange type information */
  fCommManager.AllGather(fEmbeddingStiffMessageID, fEmbeddingStiff);

  // OTHER  //
  /* reset the rhop * r array */
  frhop_r_man.SetMajorDimension(nnd, true);
  frhop_r.Resize(nnd);
  /* exchange type information */
  fCommManager.AllGather(frhop_rMessageID, frhop_r);
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
		  rp_i[k] += rhop_i * r_ij[k]; 
		  rp_j[k] += rhop_j * r_ij[k]; 
		}
	    }
	}
}


void EAMT::GetEmbEnergy(const dArray2DT& coords,const dArray2DT rho,
		  dArray2DT& Emb)
{
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
  EAMPropertyT::EmbedForceFunction emb_force = NULL;
  iArrayT neighbors;
  Emb = 0.0;
  
  for (int i = 0; i < fNeighbors.MajorDim(); i++)
    {
      fNeighbors.RowAlias(i, neighbors);
      
      int   tag_i = neighbors[0]; 
      int  type_i = fType[tag_i];

      emb_force  = fEAMProperties[type_i]->getEmbedForce();
      Emb(tag_i,0) = emb_force(rho(tag_i,0),NULL,NULL); 
    }  
}

void EAMT::GetEmbStiff(const dArray2DT& coords,const dArray2DT rho,
		       dArray2DT& Emb)
{
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
