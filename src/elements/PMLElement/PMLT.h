/* $Id: PMLT.h,v 1.9 2003-01-29 07:34:27 paklein Exp $ */
#ifndef _PML_T_H_
#define _PML_T_H_

#include "SolidElementT.h"

namespace Tahoe {

class PMLT: public SolidElementT
{
  public:
      
	PMLT(const ElementSupportT& support, const FieldT& field);

	/** initialization. called immediately after constructor */
	virtual void Initialize(void);
	/** total strain */
	/* accessors */
		
	virtual void AddNodalForce(const FieldT& field, int node, dArrayT& force);	
	virtual void AddLinearMomentum(dArrayT& momentum);
	virtual void SendOutput(int kincode);


	const dMatrixT& GradU(void) const;
	const dMatrixT& GradU(int ip) const;

	/** total strain from the end of the previous time step */
	const dMatrixT& GradU_last(void) const;
	const dMatrixT& GradU_last(int ip) const;

	/** returns 0 since total displacement is decomposed into nsd fields */
	virtual int InterpolantDOFs(void) const { return 0; };

	/** construct the field.
	 * \param nodes list of nodes at which the field should be constructed
	 * \param DOFs array of the field degrees of freedom for all nodes */
	virtual void NodalDOFs(const iArrayT& nodes, dArray2DT& DOFs) const;

  protected:
	/* construct the effective mass matrix */
	virtual void LHSDriver(GlobalT::SystemTypeT);
	void ElementLHSDriver(void);

	/* form the residual force vector */
	virtual void RHSDriver(void);
	void ElementRHSDriver(void);

	virtual void FormMass(int mass_type, double constM);
	void FormDamping(int mass_type, double constC);
	virtual void FormStiffness(double constK);

	
	virtual void FormMa(int mass_type, double constM, const LocalArrayT& body_force);
	void FormCv_PML(int mass_type, double constC, const LocalArrayT& body_force);
	virtual void FormKd(double constK);

	virtual void ReadMaterialData(ifstreamT& in);
	virtual void SetGlobalShape(void);
	/* construct output labels array */
	virtual void SetNodalOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void SetElementOutputCodes(IOBaseT::OutputModeT mode, const iArrayT& flags,
		iArrayT& counts) const;
	virtual void GenerateOutputLabels(const iArrayT& n_counts, ArrayT<StringT>& n_labels, 
		const iArrayT& e_counts, ArrayT<StringT>& e_labels) const;	
	virtual void ComputeOutput(const iArrayT& n_codes, dArray2DT& n_values,
	                           const iArrayT& e_codes, dArray2DT& e_values);

  	           
  private:
  
  	GlobalT::SystemTypeT TangentType(void) const;
	
	LocalArrayT& TotalDisp(LocalArrayT& disp);
	LocalArrayT& TotalLastDisp(LocalArrayT& disp);
	LocalArrayT& TotalVel(LocalArrayT& vel);
	LocalArrayT& TotalAcc(LocalArrayT& acc);

//Q: Who calls these accessors?
//	const LocalArrayT& LastDisplacements(void) const;
//	const LocalArrayT& Velocities(void) const;
//	const LocalArrayT& Accelerations(void) const;
	
	void Ba(dMatrixT& Ba_matrix, dMatrixT& B_matrix);
	void Bb(dMatrixT& Bb_matrix, dMatrixT& B_matrix);
	
  private:

	/** indicies of elements in the list of material needs */
	enum MaterialNeedsT {kstrain = 0,
	                kstrain_last = 1};
	int fNeedsOffset; //NOTE - better to have this or a separate array?
  
  	/** return values */
  	ArrayT<dMatrixT> fGradU_List;
  	ArrayT<dMatrixT> fGradU_last_List;
  	
  	double fdt;
  	
  	int fNEESub;
  	
  	/** work space */
  	dMatrixT fGradU;
  	dMatrixT fBa;
  	dMatrixT fBb;
  	
  	dMatrixT fDa;
  	dMatrixT fDb;

	ElementMatrixT fLHSa;
	ElementMatrixT fLHSb;
	
	dArrayT fRHSa;
	dArrayT fRHSb;
		
	LocalArrayT fTotDisp;
	LocalArrayT fTotLastDisp;	         
	LocalArrayT fTotVel;
	LocalArrayT fTotAcc;  
    
};

/* inlines */

inline const dMatrixT& PMLT::GradU(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
	{
		cout << "\n PMLT::LinearStrain: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	return fGradU_List[CurrIP()];
}

inline const dMatrixT& PMLT::GradU(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain])
	{
		cout << "\n PMLT::LinearStrain: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	return fGradU_List[ip];
}

inline const dMatrixT& PMLT::GradU_last(void) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
	{
		cout << "\n PMLT::Linearstrain last: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	return fGradU_last_List[CurrIP()];
}

inline const dMatrixT& PMLT::GradU_last(int ip) const
{
#if __option(extended_errorcheck)
	/* check need */
	int mat_num = CurrentElement().MaterialNumber();
	const ArrayT<bool>& needs = fMaterialNeeds[mat_num];
	if (!needs[fNeedsOffset + kstrain_last])
	{
		cout << "\n PMLT::LinearStrain_last: material " << mat_num + 1 
		     << " did not specify this need" << endl;
		throw ExceptionT::kGeneralFail;
	}
#endif

	return fGradU_last_List[ip];
}

} // namespace Tahoe 
#endif /* _PML_T_H_ */
