/* $Id: XDOF_ManagerT.cpp,v 1.11 2002-10-20 22:49:30 paklein Exp $ */
/* created: paklein (06/01/1998) */
/* base class which defines the interface for a manager */
/* of DOF's comprised of FE DOF's plus constrain DOF's */

#include "XDOF_ManagerT.h"

#include "DOFElementT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"

/* constructor */

using namespace Tahoe;

XDOF_ManagerT::XDOF_ManagerT(void): 
	fNumTags(0),
	fStartTag(-1) 
{ 

}

/* destructor */
XDOF_ManagerT::~XDOF_ManagerT(void)
{
	/* consistency check */
	if (fXDOF_Eqnos.Length() != fXDOFs.Length()) throw ExceptionT::kGeneralFail;

	for (int i = 0; i < fXDOF_Eqnos.Length(); i++)
	{
		delete fXDOF_Eqnos[i];
		delete fXDOFs[i];
	}
}

/* add element group to list */
void XDOF_ManagerT::XDOF_Register(DOFElementT* group, const iArrayT& numDOF)
{
	/* check start tag */
	if (fStartTag == -1)
	{
		cout << "\n XDOF_ManagerT::XDOF_Register: start tag has not been set: "
		     << fStartTag << endl;
		throw ExceptionT::kGeneralFail;
	}

	/* can only register once */
	if (!fDOFElements.AppendUnique(group)) 
	{
		cout << "\n XDOF_ManagerT::XDOF_Register: group is already registered, requesting:\n" 
		     << numDOF << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* keep number of tag sets for each group */
	fNumTagSets.Append(numDOF.Length());

	/* add to lists */
	for (int i = 0; i < numDOF.Length(); i++)
	{
		/* initialize set length */
		fTagSetLength.Append(0);

		fXDOF_Eqnos.Append(new iArray2DT(0, numDOF[i]));
		fXDOFs.Append(new dArray2DT(0, numDOF[i]));
	}
	
	/* reset number of tags */
	int group_num = fDOFElements.Length() - 1;
	if (group_num == 0)	fNumTags = fStartTag;

	/* set contact configuration */	
	ConfigureElementGroup(group_num, fNumTags);
}

/* get equation numbers for the specified constraint tags */
const iArray2DT& XDOF_ManagerT::XDOF_Eqnos(const DOFElementT* group, 
	int tag_set) const
{
	return *(fXDOF_Eqnos[TagSetIndex(group, tag_set)]);
}

/* returns reference to current constraint values */
const dArray2DT& XDOF_ManagerT::XDOF(const DOFElementT* group, int tag_set) const
{
	return *(fXDOFs[TagSetIndex(group, tag_set)]);
}

/**********************************************************************
* Protected
**********************************************************************/

/* return the number of XDOF equations in the specified group */
int XDOF_ManagerT::NumEquations(int group) const
{
	int neq = 0;
	for (int j = 0; j < fDOFElements.Length(); j++)
		if (fDOFElements[j]->Group() == group)
		{
			int num_sets = fNumTagSets[j];
			for (int k = 0; k < num_sets; k++)
			{
				int dex = TagSetIndex(fDOFElements[j], k);
				neq += fXDOF_Eqnos[dex]->Length(); /* all active */
			}
		}
	return neq;
}

/* prompt element groups to reset tags */
bool XDOF_ManagerT::ResetTags(int group)
{
	if (fDOFElements.Length() > 0)
	{
		/* query (all) groups to reconfigure */
		int relax = 0;

		/* loop over element in the group */
		int index = 0;
		for (int j = 0; j < fDOFElements.Length(); j++)
			if (fDOFElements[j]->Group() == group && 
			    fDOFElements[j]->Reconfigure() == 1)
				relax = 1;
	
		if (relax)
		{
			/* check */
			if (fStartTag == -1)
			{
				cout << "\n XDOF_ManagerT::ResetTags: start tag has not been set: "
				     << fStartTag << endl;
				throw ExceptionT::kGeneralFail;
			}
		
			/* reset */
			fNumTags = fStartTag;
		
			/* loop over DOF element groups */
			for (int i = 0; i < fDOFElements.Length(); i++)
				if (fDOFElements[i]->Group() == group)
				{
					/* reset contact configs */
					ConfigureElementGroup(i, fNumTags);

					/* restore values */
					/* loop over tag sets */
					for (int j = 0; j < fNumTagSets[i]; j++)
						fDOFElements[i]->ResetDOF(*fXDOFs[index++], j);
				}

			return true;
		}
		else
			return false;
	}
	else	
		return false;
}

/* call groups to reset external DOF's */
void XDOF_ManagerT::Reset(int group)
{
	/* loop over element */
	int index = 0;
	for (int i = 0; i < fDOFElements.Length(); i++)
	{
		if (fDOFElements[i]->Group() == group)
		{
			/* loop over tag sets */
			for (int j = 0; j < fNumTagSets[i]; j++)
				fDOFElements[i]->ResetDOF(*fXDOFs[index++], j);
		}
		else /* next element */
			index += fNumTagSets[i];
	}
}

/* update DOF's using the global update vector */
void XDOF_ManagerT::Update(int group, const dArrayT& update)
{
	/* loop over elements */
	int dex = 0;
	for (int i = 0; i < fDOFElements.Length(); i++)
	{
		int nsets = fNumTagSets[i];
		if (fDOFElements[i]->Group() == group)
		{
			for (int j = 0; j < nsets; j++)
			{
				/* group data */
				iArray2DT& XDOF_eqnos = *fXDOF_Eqnos[dex];
				dArray2DT& XDOF = *fXDOFs[dex];

				int   num_eq = XDOF_eqnos.Length();
				int*     peq = XDOF_eqnos.Pointer();
				double* pDOF = XDOF.Pointer();
		
				/* assume mapped sequentially and all are active */
				for (int k = 0; k < num_eq; k++)
					*pDOF++ += update[*peq++ - 1]; //OFFSET
			
				/* next tag set */
				dex++;
			}
		}
		else
			/* next group */
			dex += nsets;
	}
}

/* (self-)configure element group */
void XDOF_ManagerT::ConfigureElementGroup(int group_number, int& tag_num)
{
	/* first tag in set */
	int index = 0;
	for (int i = 0; i < group_number; i++)
		index += fNumTagSets[i];
		
	/* reset DOF tags arrays */	
	fDOFElements[group_number]->SetDOFTags();

	/* loop over tag sets */
	int num_sets = fNumTagSets[group_number];
	for (int set = 0; set < num_sets; set++)
	{
		/* get current DOF tags array */
		iArrayT& DOF_tags = fDOFElements[group_number]->DOFTags(set);
		fTagSetLength[set] = DOF_tags.Length();
							
		/* resize global arrays */
		fXDOF_Eqnos[index]->Resize(DOF_tags.Length(), false);
		fXDOFs[index]->Resize(DOF_tags.Length(), false);

		/* initialize */
		*(fXDOFs[index]) = 0.0;
			
		/* generate new DOF tags */
		for (int k = 0; k < DOF_tags.Length(); k++)
			DOF_tags[k] = tag_num++;
	
		/* next set */
		index++;
	}
	
	/* call to (self-)configure */
	fDOFElements[group_number]->GenerateElementData();
}

/* assign equation numbers */
void XDOF_ManagerT::SetEquations(int group, int& num_eq)
{
	/* loop over elements */
	int dex = 0;
	for (int i = 0; i < fDOFElements.Length(); i++)
	{
		int nsets = fNumTagSets[i];
		if (fDOFElements[i]->Group() == group)
		{
			/* assign number sequentially through the set */
			for (int j = 0; j < nsets; j++)
			{
				int* peqno = fXDOF_Eqnos[dex]->Pointer();
				int length = fXDOF_Eqnos[dex]->Length();
				for (int k = 0; k < length; k++)
					*peqno++ = ++num_eq;

				/* next tag set */
				dex++;
			}
		}
		else
			/* next group */
			dex += nsets;
	}
}

/* access to global equation numbers */
void XDOF_ManagerT::EquationNumbers(int group, AutoArrayT<iArray2DT*>& equationsets)
{
	/* loop over elements */
	for (int i = 0; i < fDOFElements.Length(); i++)
		if (fDOFElements[i]->Group() == group)
			equationsets.Append(fXDOF_Eqnos[i]);
}

/* remove external DOF's from first slot of each row */
void XDOF_ManagerT::CheckEquationNumbers(ostream& out, iArray2DT& eqnos)
{
//NOTE: The augmented Lagrangian formulation puts a zero on the
//      diagonal of the unfactorized matrix. This might be OK if
//      that equation isn't exactly the first one. Check this and
//      swap equations with any of the displacement DOF's connected
//      to the aug. lag. DOF.

	/* quick exit */
	if (fDOFElements.Length() == 0) return;

	out << "\n XDOF_ManagerT::CheckEquationNumbers: swapped node/tag's:\n\n";

	/* output formatting */
	const int num_cols = 3;
	int col = 0;

	/* put external equation number after displacement equations */
	int set_index =-1;
	for (int i = 0; i < fDOFElements.Length(); i++)
		for (int set = 0; set < fNumTagSets[i]; set++)
		{
			/* next */
			set_index++;
			
			iArray2DT& XDOF_eqnos = *fXDOF_Eqnos[set_index];
			if (XDOF_eqnos.MinorDim() == 1)
			{
				iArrayT eq_temp;
				eq_temp.Alias(XDOF_eqnos);

				/* constraint connectivities */
				const iArray2DT& connects = fDOFElements[i]->DOFConnects(set);
				iArrayT elem;
				iArrayT eq_u;
				int nen = connects.MinorDim();
				int x_node = 0;
				for (int j = 0; j < connects.MajorDim(); j++)
				{
					connects.RowAlias(j, elem);
					if (elem[nen-1] < 0 ) 
					{ 	/* no xdof exists for this u node */
						out << "\n XDOF_ManagerT: found u node w/o xdof: ";
                        out << elem[nen-1] << endl;

					}
					else if (elem[nen-1] + 1 <= eqnos.MajorDim()) 
					{	/* xdof node is the u node range */
						cout << "\n XDOF_ManagerT: expecting external DOF tag in final slot: ";
						cout << elem[nen-1] << endl;
						throw ExceptionT::kGeneralFail;
					} else {
		
						/* get equations from 1st node in connect */
						eqnos.RowAlias(elem[0], eq_u); // connects are 1,...

						/* find highest equation of the node */
						int max_eq = eq_u.Max();
						if (max_eq < 0)
						{	/* xdof has been assigned to a node with ess.bcs*/
							out << "\n XDOF_ManagerT: no active u equations found with xdof tag ";
							out << elem[nen-1] << '\n';
							out << " connectivity: ";
							out << elem.no_wrap() << '\n';
							out << endl;
						} else {

							/* swap */
							eq_u.ChangeValue(max_eq, XDOF_eqnos(x_node,0));
							XDOF_eqnos(x_node,0) = max_eq;
		
							out << setw(kIntWidth) << elem[0];
							out << setw(kIntWidth) << elem[nen-1];
							if (++col == num_cols)
							{
								out << '\n';
								col = 0;
							}
							else
							out << setw(kIntWidth) << " ";
	
						}
						x_node++;
					}
				}
			}
			else
			{
				cout << "\n XDOF_ManagerT::CheckNewEquationNumbers: no check for NDOF > 1: " 
				     << XDOF_eqnos.MinorDim() << endl;
			}
		}
	
	/* mid-line */
	if (col != 0) out << '\n';
}

/* resolve index of the tag set */
int XDOF_ManagerT::TagSetIndex(const DOFElementT* group, int tag_set) const
{
	/* find matching element group */
	int index = fDOFElements.PositionOf((DOFElementT*) group);
	if (index < 0)
	{
		cout << "\n XDOF_ManagerT::TagSetIndex: no matching group found" << endl;
		throw ExceptionT::kGeneralFail;
	}
	
	/* check tag set */
	if (tag_set < 0 || tag_set >= fNumTagSets[index])
	{
		cout << "\n XDOF_ManagerT::TagSetIndex: tag set number " << tag_set 
		     << " is out of range {0, " << fNumTagSets[index] - 1 << "}" << endl;
		throw ExceptionT::kOutOfRange;
	}
	
	/* offset to groups tag sets */
	int offset = 0;
	for (int i = 0; i < index; i++)
		offset += fNumTagSets[index];

	/* resolve set index */
	offset += tag_set;
	if (offset >= fXDOF_Eqnos.Length())
	{
		cout << "\n XDOF_ManagerT::TagSetIndex: error resolving tag set index: " 
		     << offset << " > " << fXDOF_Eqnos.Length() << endl;
		throw ExceptionT::kGeneralFail;
	}	
	return offset;
}

/* resolve tag into its tag set and tag offset */
bool XDOF_ManagerT::ResolveTagSet(int tag, int& tag_set, int& tag_set_start) const
{
	tag_set = tag_set_start = -1;
	int offset = fStartTag;
	bool found = false;
	for (int i = 0; i < fTagSetLength.Length() && !found; i++)
	{
		/* in range */
		if (tag - offset < fTagSetLength[i])
		{
			found = true;
			tag_set = i;
			tag_set_start = offset;
		}
		else /* next set - assuming tags assigned to sets sequentially */
			offset += fTagSetLength[i];
	}
	return found;
}
