/* $Id: XDOF_ManagerT.cpp,v 1.4 2001-08-29 07:08:22 paklein Exp $ */
/* created: paklein (06/01/1998) */
/* base class which defines the interface for a manager */
/* of DOF's comprised of FE DOF's plus constrain DOF's */

#include "XDOF_ManagerT.h"

#include "DOFElementT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"

/* constructor */
XDOF_ManagerT::XDOF_ManagerT(void): 
	fNumTags(0),
	fStartTag(-1) 
{ 

}

/* destructor */
XDOF_ManagerT::~XDOF_ManagerT(void)
{
	/* consistency check */
	if (fXDOF_Eqnos.Length() != fXDOFs.Length()) throw eGeneralFail;

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
		throw eGeneralFail;
	}

	/* can only register once */
	if (!fDOFElements.AppendUnique(group)) 
	{
		cout << "\n XDOF_ManagerT::XDOF_Register: group is already registered, requesting:\n" 
		     << numDOF << endl;
		throw eGeneralFail;
	}
	
	/* keep number of tag sets for each group */
	fNumTagSets.Append(numDOF.Length());

	/* initialize set length */
	fTagSetLength.Append(0);

	/* add to lists */
	for (int i = 0; i < numDOF.Length(); i++)
	{
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

/* prompt element groups to reset tags */
bool XDOF_ManagerT::ResetTags(void)
{
	if (fDOFElements.Length() > 0)
	{
		/* query (all) groups to reconfigure */
		int relax = 0;
		for (int j = 0; j < fDOFElements.Length(); j++)
		{
			if (fDOFElements[j]->Reconfigure() == 1)
				relax = 1;
		}
	
		if (relax)
		{
			/* check */
			if (fStartTag == -1)
			{
				cout << "\n XDOF_ManagerT::ResetTags: start tag has not been set: "
				     << fStartTag << endl;
				throw eGeneralFail;
			}
		
			/* reset */
			fNumTags = fStartTag;
		
			/* loop over DOF element groups */
			for (int i = 0; i < fDOFElements.Length(); i++)
			{
				/* reset contact configs */
				ConfigureElementGroup(i, fNumTags);
				
				/* restore values */
				fDOFElements[i]->ResetDOF(*fXDOFs[i], 0);
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
void XDOF_ManagerT::Reset(void)
{
	/* loop over groups */
	int index = 0;
	for (int i = 0; i < fDOFElements.Length(); i++)
		/* loop over tag sets */
		for (int j = 0; j < fNumTagSets[i]; j++)
			fDOFElements[i]->ResetDOF(*fXDOFs[index++], j);
}

/* update DOF's using the global update vector */
void XDOF_ManagerT::Update(const dArrayT& update)
{
	/* loop over all tag sets */
	for (int i = 0; i < fXDOF_Eqnos.Length(); i++)
	{
		/* group data */
		iArray2DT& XDOF_eqnos = *fXDOF_Eqnos[i];
		dArray2DT& XDOF = *fXDOFs[i];

		int   num_eq = XDOF_eqnos.Length();
		int*     peq = XDOF_eqnos.Pointer();
		double* pDOF = XDOF.Pointer();
		
		/* assume mapped sequentially and all are active */
		for (int j = 0; j < num_eq; j++)
			*pDOF++ += update[*peq++ - 1]; //OFFSET
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
void XDOF_ManagerT::SetEquations(int& num_eq)
{
	int num_sets = fXDOF_Eqnos.Length();
	for (int i = 0; i < num_sets; i++)
	{
		int* peqno = fXDOF_Eqnos[i]->Pointer();
		int length = fXDOF_Eqnos[i]->Length();
		for (int j = 0; j < length; j++)
			*peqno++ = ++num_eq;
	}
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
			if (XDOF_eqnos.MinorDim() > 1)
			{
				cout << "\n XDOF_ManagerT::CheckNewEquationNumbers: expecting only 1 DOF" << endl;
				throw eGeneralFail;
			}
			iArrayT eq_temp;
			eq_temp.Alias(XDOF_eqnos);

			/* constraint connectivities */
			const iArray2DT& connects = fDOFElements[i]->DOFConnects(set);
			iArrayT elem;
			iArrayT eq_u;
			int nen = connects.MinorDim();
			for (int j = 0; j < connects.MajorDim(); j++)
			{
				connects.RowAlias(j, elem);
				if (elem[nen-1] + 1 <= eqnos.MajorDim()) /* eqnos has disp DOF's */
				{
					cout << "\n XDOF_ManagerT: expecting external DOF tag in final slot: ";
					cout << elem[nen-1] << endl;
					throw eGeneralFail;
				}
		
				/* get equations from 1st node in connect */
				eqnos.RowAlias(elem[0], eq_u); // connects are 1,...

				/* find highest equation of the node */
				int max_eq = eq_u.Max();
				if (max_eq < 0)
				{
					cout << "\n XDOF_ManagerT: no active equations found with contact tag ";
					cout << elem[nen-1] << '\n';
					cout << " connectivity: ";
					cout << elem.no_wrap() << '\n';
					cout << endl;
					throw eGeneralFail;
				}
			
				/* swap */
				eq_u.ChangeValue(max_eq, XDOF_eqnos(j,0));
				XDOF_eqnos(j,0) = max_eq;
			
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
		throw eGeneralFail;
	}
	
	/* check tag set */
	if (tag_set < 0 || tag_set >= fNumTagSets[index])
	{
		cout << "\n XDOF_ManagerT::TagSetIndex: tag set number " << tag_set 
		     << " is out of range {0, " << fNumTagSets[index] - 1 << "}" << endl;
		throw eOutOfRange;
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
		throw eGeneralFail;
	}	
	return offset;
}

/* resolve tag into its tag set and tag offset */
bool XDOF_ManagerT::ResolveTagSet(int tag, int& tag_set, int& tag_set_start)
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
