/* $Id: XDOF_ManagerT.cpp,v 1.1.1.1 2001-01-29 08:20:40 paklein Exp $ */
/* created: paklein (06/01/1998)                                          */
/* base class which defines the interface for a manager                   */
/* of DOF's comprised of FE DOF's plus constrain DOF's                    */

#include "XDOF_ManagerT.h"

#include "DOFElementT.h"
#include "iArrayT.h"
#include "iArray2DT.h"
#include "dArray2DT.h"

/* constructor */
XDOF_ManagerT::XDOF_ManagerT(void) { }

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
void XDOF_ManagerT::Register(DOFElementT* group, int numDOF)
{
	/* can only register once */
	if ( !fDOFElements.AppendUnique(group) ) throw eGeneralFail;

	/* add to lists */
	fXDOF_Eqnos.Append(new iArray2DT(0, numDOF));
	fXDOFs.Append(new dArray2DT(0, numDOF));
}

/* get equation numbers for the specified constraint tags */
const iArray2DT& XDOF_ManagerT::XDOF_Eqnos(const DOFElementT* group) const
{
	/* find matching element group */
	for (int i = 0; i < fDOFElements.Length(); i++)
		if (fDOFElements[i] == group) return *(fXDOF_Eqnos[i]);

	/* should never fall through */
	throw eGeneralFail;
	return *(fXDOF_Eqnos[0]);
}

/* returns reference to current constraint values */
const dArray2DT& XDOF_ManagerT::XDOF(const DOFElementT* group) const
{
	/* find matching element group */
	for (int i = 0; i < fDOFElements.Length(); i++)
		if (fDOFElements[i] == group) return *(fXDOFs[i]);

	/* should never fall through */
	throw eGeneralFail;
	return *(fXDOFs[0]);
}

/**********************************************************************
* Protected
**********************************************************************/

/* call groups to reset external DOF's */
void XDOF_ManagerT::Reset(void)
{
	/* reset */
	for (int i = 0; i < fDOFElements.Length(); i++)
		fDOFElements[i]->ResetDOF(*fXDOFs[i]);
}

/* update DOF's using the global update vector */
void XDOF_ManagerT::Update(const dArrayT& update)
{
	for (int i = 0; i < fDOFElements.Length(); i++)
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
	/* get current DOF tags array */
	iArrayT& DOFtags = fDOFElements[group_number]->SetDOFTags();
							
	/* resize global arrays */
	fXDOF_Eqnos[group_number]->Resize(DOFtags.Length(), false);
	fXDOFs[group_number]->Resize(DOFtags.Length(), false);

	/* initialize */
	*(fXDOFs[group_number]) = 0.0;
			
	/* generate new DOF tags */
	for (int k = 0; k < DOFtags.Length(); k++)
		DOFtags[k] = tag_num++;
	
	/* call to (self-)configure */
	fDOFElements[group_number]->GenerateElementData();
}

/* assign equation numbers */
void XDOF_ManagerT::SetEquations(int& num_eq)
{
	int num_groups = fXDOF_Eqnos.Length();
	for (int i = 0; i < num_groups; i++)
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

	out << "\n XDOF_ManagerT::CheckEquationNumbers: swapped node/tag's:\n\n";

	/* output formatting */
	const int num_cols = 3;
	int col = 0;

	/* put external equation number after displacement equations */
	for (int i = 0; i < fDOFElements.Length(); i++)
	{
		iArray2DT& XDOF_eqnos = *fXDOF_Eqnos[i];
		if (XDOF_eqnos.MinorDim() > 1)
		{
			cout << "\n XDOF_ManagerT::CheckNewEquationNumbers: expecting only 1 DOF" << endl;
			throw eGeneralFail;
		}
		iArrayT eq_temp;
		eq_temp.Alias(XDOF_eqnos);

		/* constraint connectivities */
		const iArray2DT& connects = fDOFElements[i]->DOFConnects();
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
