/* $Id: XDOF_ManagerT.h,v 1.9 2004-01-05 07:12:48 paklein Exp $ */
/* created: paklein (06/01/1998) */
#ifndef _XDOF_MANAGER_T_H_
#define _XDOF_MANAGER_T_H_

/* direct members */
#include "dArrayT.h"
#include "AutoArrayT.h"
#include "VariArrayT.h"
#include "GlobalT.h"

namespace Tahoe {

/* forward declarations */
class DOFElementT;
class iArray2DT;
class dArray2DT;
class iArrayT;
template <class TYPE> class RaggedArray2DT;

/** mix-in class for manager of degrees of freedom requested
 * by DOFElementT's. Element groups must be derived from the
 * DOFElementT class and must register themselves using the
 * XDOF_ManagerT::Register funciton. */
class XDOF_ManagerT
{
public:
	
	/** constructor */
	XDOF_ManagerT(void);	

	/** destructor */
	virtual ~XDOF_ManagerT(void);

	/** add element group to list. Each group can request an arbitrary
	 * number of tag sets. Each set can have a different number of DOF's
	 * per tag requested. Requests for tags are collected by the XDOF_ManagerT
	 * using the DOFElementT interface.
	 * \param group pointer to the DOFElementT class. Each group should
	 *        only register once.
	 * \param numDOF array of the number of degrees of freedom per tag 
	 *        in each set of tags. The length of the array is the number
	 *        of tag sets the group requires. */
	virtual void XDOF_Register(DOFElementT* group, const iArrayT& numDOF);
	
	/** get equations of the element DOF's.
	 * \param group pointer to the DOFElementT requesting equation numbers
	 * \param tag_set set number with the DOFElementT */
	virtual const iArray2DT& XDOF_Eqnos(const DOFElementT* group, int tag_set) const;

	/** get values of the element DOF's.
	 * \param group pointer to the DOFElementT requesting DOF values
	 * \param tag_set set number with the DOFElementT */
	virtual const dArray2DT& XDOF(const DOFElementT* group, int tag_set) const;

	/** \name collecting equation numbers 
	 * These methods foresee the need to collect equation numbers using connectivities
	 * that contain both XDOF tags and node numbers */
	/*@{*/
	/** collection equation numbers for mixed connectivity. Connectivity
	 * can be either node numbers of tag numbers obtained through the
	 * inherited XDOF_ManagerT interface. For tags that are nodes, all
	 * equations for that node across all fields in the group.
	 * \param group equation group number
	 * \param nodes element connectivity: [nen]
	 * \param eqnos destination for equation numbers: [nen] x [ndof_i] */
	virtual void XDOF_SetLocalEqnos(int group, const iArrayT& nodes, iArray2DT& eqnos) = 0;

	/** collection equation numbers for mixed connectivities. Connectivities
	 * can be either node numbers of tag numbers obtained through the
	 * inherited XDOF_ManagerT interface. For tags that are nodes, all
	 * equations for that node across all fields in the group.
	 * \param group equation group number
	 * \param nodes element connectivities: [nel] x [nen]
	 * \param eqnos destination for equation numbers: [nel] x [nen*ndof_j] */
	virtual void XDOF_SetLocalEqnos(int group, const iArray2DT& nodes, iArray2DT& eqnos) const = 0;

	/** collection equation numbers for mixed connectivities. Connectivities
	 * can be either node numbers of tag numbers obtained through the
	 * inherited XDOF_ManagerT interface. Connectivities are passed in
	 * a RaggedArray2DT, which allows an arbitrary number of nodes per
	 * element. For tags that are nodes, all
	 * equations for that node across all fields in the group.
	 * \param group equation group number
	 * \param nodes element connectivities: [nel] x [nen_i]
	 * \param eqnos destination for equation numbers: [nel] x [nen_i*ndof_j] */
	virtual void XDOF_SetLocalEqnos(int group, const RaggedArray2DT<int>& nodes, 
		RaggedArray2DT<int>& eqnos) const = 0;
	/*@}*/

protected:

	/** return the number of XDOF equations in the specified group */
	int NumEquations(int group) const;

	/** set the start tag */
	void SetStartTag(int start_tag) { fStartTag = start_tag; };

	/** prompt elements to restore their internal states */
	void ResetState(int group);

	/** prompt elements in the specified group to reset tags.
	 * \return true if tags have been reset */
	GlobalT::RelaxCodeT ResetTags(int group);

	/** return the total number of tag sets */
	int NumTagSets(void) const { return fXDOF_Eqnos.Length(); };

	/** call elements in specified group the to reset external DOF's */
	void Reset(int group);

	/** update DOF's in the specified group using the global 
	 * update vector */
	void Update(int group, const dArrayT& update);

	/** (self-)configure element group */
	void ConfigureElementGroup(int group_number, int& tag_num);

	/** assign equation numbers */
	void SetEquations(int group, int& num_eq);

	/** remove external DOF's from first slot of each row. The augmented Lagrangian 
	 * formulation puts a zero on the diagonal of the unfactorized matrix. This might 
	 * be OK if that equation isn't exactly the first one. Check this and swap equations 
	 * with any of the displacement DOF's connected to the augmented Lagrangian DOF. */
	void CheckEquationNumbers(ostream& out, iArray2DT& eqnos);

	/** append equation numbers for the tags sets in the specified group */
	void EquationNumbers(int group, AutoArrayT<iArray2DT*>& equationsets);
	
	/** resolve index of the tag set.
	 * \param group pointer to the DOFElementT
	 * \param tag_set set number for the element 
	 * \return the index of the tag set in XDOF_ManagerT::fXDOF_Eqnos
	 *        and XDOF_ManagerT::fXDOF */
	int TagSetIndex(const DOFElementT* group, int tag_set) const;

	/** resolve tag into its tag set and tag offset.
	 * \param tag tag to resolve
	 * \param tag_set returns with the set containing the specified tag
	 * \param tag_set_start first tag number in the set
	 * \return true of the tag is successfully resolved, false otherwise. 
	 * \note The alorithm to determine the tag set assumes the tags are
	 *       assigned in the order the sets appear in XDOF_ManagerT::fXDOF_Eqnos
	 *       and XDOF_ManagerT::fXDOF. */
	bool ResolveTagSet(int tag, int& tag_set, int& tag_set_start) const;

protected:

	/** registered element groups */
	AutoArrayT<DOFElementT*> fDOFElements;
	
	/** \name tag info */
	/*@{*/
	int fStartTag; /**< first tag number */
	int fNumTags;  /**< total number of tags */
	AutoArrayT<int> fNumTagSets;   /**< number of tag sets for each element group */
	AutoArrayT<int> fTagSetLength; /**< number of tags in each tag set */
	/*@}*/

	/** global equations numbers for each tag set */
	AutoArrayT<iArray2DT*> fXDOF_Eqnos;

	/** values of the degrees of freedom by set */
	AutoArrayT<dArray2DT*> fXDOFs;
};

} // namespace Tahoe 
#endif /* _XDOF_MANAGER_T_H_ */
