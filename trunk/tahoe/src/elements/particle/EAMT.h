/* $Id: EAMT.h,v 1.8 2003-05-08 00:44:19 saubry Exp $ */
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
	EAMT(const ElementSupportT& support, const FieldT& field);

	/** collecting element group equation numbers */
	virtual void Equations(AutoArrayT<const iArray2DT*>& eq_1,
		AutoArrayT<const RaggedArray2DT<int>*>& eq_2);

	/** class initialization. Among other things, element work space
	 * is allocated and connectivities are read. */
	virtual void Initialize(void);

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

protected:

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

	/** construct the list of properties from the given input stream */
	virtual void EchoProperties(ifstreamT& in, ofstreamT& out);

	/** generate labels for output data */
	virtual void GenerateOutputLabels(ArrayT<StringT>& labels) const;

private:

	void GetRho2D(const dArray2DT& coords,dArray2DT& rho);
	void GetRho3D(const dArray2DT& coords,dArray2DT& rho);

	void GetRhoForce2D(const dArray2DT& coords,dArray2DT& rho);
	void GetRhoForce3D(const dArray2DT& coords,dArray2DT& rho);

	void GetRhoStiff2D(const dArray2DT& coords,dArray2DT& rho);
	void GetRhoStiff3D(const dArray2DT& coords,dArray2DT& rho);

	void GetEmbEnergy(const dArray2DT& coords,const dArray2DT rho,
			  dArray2DT& Emb);
	void GetEmbForce(const dArray2DT& coords,const dArray2DT rho,
			 dArray2DT& Emb);

	void EAMT::GetEmbStiff(const dArray2DT& coords,const dArray2DT rho,
			       dArray2DT& Emb);

	/** particle pair-properties list */
	ArrayT<EAMPropertyT*> fEAMProperties;

	/** neighbor lists */
	RaggedArray2DT<int> fNeighbors;

	/** equation numbers */
	RaggedArray2DT<int> fEqnos;

	/** electron density */
	dArray2DT fElectronDensity;
	nVariArray2DT<double> fElectronDensity_man;
	int fElectronDensityMessageID;

	/** electron density force*/
	dArray2DT fElectronDensityForce;
	nVariArray2DT<double> fElectronDensityForce_man;
	int fElectronDensityForceMessageID;

	/** electron density stiffness*/
	dArray2DT fElectronDensityStiff;
	nVariArray2DT<double> fElectronDensityStiff_man;
	int fElectronDensityStiffMessageID;

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


	/** \name workspace for EAMT::RHSDriver. Used to accumulate the force for
	 * a single row of EAMT::fNeighbors. */
	/*@{*/
	dArrayT fForce_list;
	VariArrayT<double> fForce_list_man;

	/** constant matrix needed to compute the stiffness */
	dMatrixT fOneOne;
	/*@}*/
};

} /* namespace Tahoe */

#endif /* _EAM_T_H_ */

