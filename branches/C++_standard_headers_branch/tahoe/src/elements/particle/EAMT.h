/* $Id: EAMT.h,v 1.24 2005-04-08 16:41:48 d-farrell2 Exp $ */

#ifndef _EAM_T_H_
#define _EAM_T_H_

/* base class */
#include "ParticleT.h"

/* direct members */
#include "RaggedArray2DT.h"
#include "VariArrayT.h"

namespace Tahoe {

/* forward declarations */
class EAMPropertyT;

/** base class for particle types */
class EAMT: public ParticleT
{
public:

	/** constructor */
	EAMT(const ElementSupportT& support);

	/** collecting element group equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** \name connectivities.
	 * See ElementBaseT::ConnectsX and ElementBaseT::ConnectsU for more
	 * information about what these are used for */
	/*@{*/
	/** collecting element geometry connectivities */
	virtual void ConnectsX(AutoArrayT<const iArray2DT*>& connects) const;

	/** collecting element field connectivities */
	virtual void ConnectsU(AutoArrayT<const iArray2DT*>& connects_1,
	             AutoArrayT<const RaggedArray2DT<int>*>& connects_2) const;
	/*@}*/

	/** write output. Writes output data for all tags listed in
	 * ParticleT::fGlobalTag. The values per node are those specified
	 * by EAMT::GenerateOutputLabels. */
	virtual void WriteOutput(void);

	/** compute the part of the stiffness matrix associated with rows of the given 
	 * particles. This is the mixed part of the stiffness matrix involving free
	 * particles and ghost particles which have prescribed motion. */
	virtual void FormStiffness(const InverseMapT& col_to_col_eq_row_map,
		const iArray2DT& col_eq, dSPMatrixT& stiffness);

	/** set external electron density pointers */
	void SetExternalElecDensity(const dArray2DT& elecdens, const iArrayT& ghostatoms);
	
	/** set external embedding force pointers */
	void SetExternalEmbedForce(const dArray2DT& embforce, const iArrayT& ghostatoms);

	/** \name implementation of the ParameterInterfaceT interface */
	/*@{*/

	/** information about subordinate parameter lists */
	virtual void DefineSubs(SubListT& sub_list) const;

	/** return the description of the given inline subordinate parameter list */
	virtual void DefineInlineSub(const StringT& name, ParameterListT::ListOrderT& order, 
		SubListT& sub_lists) const;

	/** a pointer to the ParameterInterfaceT of the given subordinate */
	virtual ParameterInterfaceT* NewSub(const StringT& name) const;

	/** accept parameter list */
	virtual void TakeParameterList(const ParameterListT& list);
	/*@}*/

protected:

	/** list/index of output values */
	enum EAMOutputCodeT {
		kDisplacement = 0,
		kPE = 1,
		kKE = 2,
		kStress = 3,
		kStrain = 4,
		kSlipVector = 5,
		kCS = 6,
		kCN = 7
	};

	/** \name drivers called by ElementBaseT::FormRHS and ElementBaseT::FormLHS */
	/*@{*/
	/** form group contribution to the LHS matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);

	/** form group contribution to the residual */
	virtual void RHSDriver(void);
	void RHSDriver2D(void);
	void RHSDriver3D(void);
	/*@}*/
	
	/** set neighborlists and any other system configuration information
	 * based on the current information. Uses ParticleT::GenerateNeighborList
	 * to determine the neighborlists. */
	virtual void SetConfiguration(void);

	/** extract the properties information from the parameter list. See ParticleT::ExtractProperties */
	virtual void ExtractProperties(const ParameterListT& list, const ArrayT<StringT>& type_names,
		ArrayT<ParticlePropertyT*>& properties, nMatrixT<int>& properties_map);

	/** return number of values for each output variable */
	void SetOutputCount(const iArrayT& flags, iArrayT& counts) const;

	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const;

	/** return a new EAM property or NULL if the name is invalid */
	EAMPropertyT* New_EAMProperty(const StringT& name, bool throw_on_fail) const;

private:

	void GetRho2D(const dArray2DT& coords,dArray2DT& rho);
	void GetRho3D(const dArray2DT& coords,dArray2DT& rho);
	void GetRhop_r(const dArray2DT& coords,dArray2DT& rho);

	void GetEmbEnergy(const dArray2DT& coords,const dArray2DT rho,
			  dArray2DT& Emb);
	void GetEmbForce(const dArray2DT& coords,const dArray2DT rho,
			 dArray2DT& Emb);
	void GetEmbStiff(const dArray2DT& coords,const dArray2DT rho,
			       dArray2DT& Emb);

private:

	/** particle pair-properties list */
	ArrayT<EAMPropertyT*> fEAMProperties;

	/** neighbor lists */
	RaggedArray2DT<int> fNeighbors;

	/** equation numbers */
	RaggedArray2DT<int> fEqnos;

	/** \name nearest neighbor lists needed for calculation slip vector
	 * and strain */
	/*@{*/
	RaggedArray2DT<int> fNearestNeighbors;
	RaggedArray2DT<int> fRefNearestNeighbors;
	/*@}*/

	/** electron density */
	dArray2DT fElectronDensity;
	nVariArray2DT<double> fElectronDensity_man;
	int fElectronDensityMessageID;

	/** embedding energy */
	dArray2DT fEmbeddingEnergy;
	nVariArray2DT<double> fEmbeddingEnergy_man;
	int fEmbeddingEnergyMessageID;

	/** embedding force */
	dArray2DT fEmbeddingForce;
	nVariArray2DT<double> fEmbeddingForce_man;
	int fEmbeddingForceMessageID;

	/** embedding stiffness */
	dArray2DT fEmbeddingStiff;
	nVariArray2DT<double> fEmbeddingStiff_man;
	int fEmbeddingStiffMessageID;

	/** rhop * r */
	dArray2DT frhop_r;
	nVariArray2DT<double> frhop_r_man;
	int frhop_rMessageID;

	/** external electron density, embedding force */
	const dArray2DT* fExternalElecDensity;
	const dArray2DT* fExternalEmbedForce;
	const iArrayT* fExternalElecDensityNodes;
	const iArrayT* fExternalEmbedForceNodes;

	/** \name workspace for EAMT::RHSDriver. Used to accumulate the force for
	 * a single row of EAMT::fNeighbors. */
	/*@{*/
	dArrayT fForce_list;
	VariArrayT<double> fForce_list_man;

	/** constant matrix needed to compute the stiffness */
	dMatrixT fOneOne;
	/*@}*/

	/** output flags */
	iArrayT fOutputFlags;
};

} /* namespace Tahoe */

#endif /* _EAM_T_H_ */

