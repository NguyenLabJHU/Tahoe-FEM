/* $Id: DOFElementT.h,v 1.3 2004-01-05 07:35:36 paklein Exp $ */
/* created: paklein (06/01/1998) */

#ifndef _DOF_ELEMENT_T_H_
#define _DOF_ELEMENT_T_H_

namespace Tahoe {

/* forward declarations */
class iArrayT;
class dArray2DT;
class iArray2DT;

/** Mix-in interface for element types that generate degrees
 * of freedom that are solved as part of the global system.
 * This interface is used by the XDOF_ManagerT class to manage
 * additional degrees of freedom that are generated at the element
 * level. To initialize generation of element degrees of freedom,
 * the element must register with the XDOF_ManagerT using
 * XDOF_ManagerT::Register. After registration, the XDOF_ManagerT
 * uses this interface to allocate/re-allocate element degrees
 * of freedom for each load increment. The number of tag sets
 * requested by the group is set during with XDOF_ManagerT::Register.
 * Initialization of the element tags is by the following sequence:\n
 * (1) DOFElementT::SetDOFTags - determine number of tags needed in each set\n
 * (2) DOFElementT::GenerateElementData - reset connectivities with new tag sets\n
 * (3) DOFElementT::ResetDOF - initialize degrees of values\n*/
class DOFElementT
{
public:

	/** constructor */
	DOFElementT(void);
	
	/** return the equation group to which the generate degrees of
	 * freedom belong. */
	virtual int Group(void) const = 0;
	
	/** determine number of tags needed. A tag is analogous to a node
	 * number. Each tag may have an arbitrary number of degrees of freedom
	 * set with the call to XDOF_ManagerT::Register. During this call, the
	 * element determine the number of tags needed and allocate a iArrayT's
	 * into which the XDOF_ManagerT will write the tags allocate to the
	 * element. After this call, the arrays returned by DOFElementT::DOFTags
	 * should be current. */
	virtual void SetDOFTags(void) = 0;
	
	/** return the array tag numbers in the specified set currently 
	 * used/need by the group.
	 * \param tag_set tag set number within the group
	 * \return array allocated to the current number of needed/used tags */
	virtual iArrayT& DOFTags(int tag_set) = 0;

	/** generate nodal connectivities. This is an indication to the group
	 * that the array returned to the XDOF_ManagerT by DOFElementT::SetDOFTags 
	 * contains a "current" list of tags. Since the sequence of setting global 
	 * equation numbers is controlled externally, responsibility for calling the 
	 * element group to (re-)configure is also left to calls from the outside to
	 * insure the tags are current. */ 
	virtual void GenerateElementData(void) = 0;

	/** return the connectivities associated with the element generated
	 * degrees of freedom. The connectivities are used by the XDOF_ManagerT
	 * to re-order the node/tag numbers in the connectivities such that
	 * the tag associated with the element degrees of freedom lies last.
	 * For degrees of freedom like Lagrange multipliers, this sequence
	 * may be required by matrix solvers that do not allow pivoting during
	 * factorization. This call may also return an empty array to skip
	 * any checking of the node/tag number sequence. */
	virtual const iArray2DT& DOFConnects(int tag_set) const = 0;

	/** restore/initialize the values of the element degrees of freedom to the 
	 * last converged solution. The major dimension of DOF will be the number
	 * of tags requested by the element for the current time increment. Note
	 * that is the element generates a changing number of tags, this major
	 * dimension will generally by different from its dimension during the
	 * increment during which the element degrees of freedom where calculated.
	 * The element is responsible for mapping the values from the previous
	 * configuration into the into their new positions in DOF. */
	virtual void ResetDOF(dArray2DT& DOF, int tag_set) const = 0;

	/** check element group for tag reconfiguration. The element should
	 * determine whether the number of tags needed for the next time
	 * increment, given the current configuration, is different from the
	 * current number of tags. Returning 1 here will results in subsequent
	 * calls to DOFElementT::SetDOFTags, DOFElementT::GenerateElementData,
	 * and DOFElementT::ResetDOF.
	 * \return 1 if group needs to reconfigure its tags, else 0 */
	virtual int Reconfigure(void) = 0;
	
	/** restore any state data to the previous converged state. This is an indication
	 * to the group that the solution for the current time step has failed. This
	 * call will be followed by a call to DOFElementT::Reconfigure so the group
	 * can determine the number of tags needed */
	virtual void ResetState(void) = 0;
};

} // namespace Tahoe 
#endif /* _DOF_ELEMENT_T_H_ */
