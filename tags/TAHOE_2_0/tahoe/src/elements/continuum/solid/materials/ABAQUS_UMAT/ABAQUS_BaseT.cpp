/* $Id: ABAQUS_BaseT.cpp,v 1.6 2004-07-15 08:26:33 paklein Exp $ */
#include "ABAQUS_BaseT.h"

#ifdef __F2C__

#include "ElementCardT.h"
#include "dMatrixT.h"
#include "dSymMatrixT.h"


#include <ctype.h>

using namespace Tahoe;

/* constructor */
ABAQUS_BaseT::ABAQUS_BaseT(void)
{

}

/***********************************************************************
 * Protected
 ***********************************************************************/

/* conversion functions */
void ABAQUS_BaseT::dMatrixT_to_ABAQUS(const dMatrixT& A,
	nMatrixT<doublereal>& B) const
{
#if __option(extended_errorcheck)
	/* expecting ABAQUS matrix to be 3D always */
	if (B.Rows() != 3 ||
	    B.Cols() != 3 ||
	    A.Rows() != A.Cols()) ExceptionT::GeneralFail("ABAQUS_BaseT::dMatrixT_to_ABAQUS");
#endif

	if (A.Rows() == 2)
	{
		B(0,0) = doublereal(A(0,0));
		B(1,0) = doublereal(A(1,0));
		B(0,1) = doublereal(A(0,1));
		B(1,1) = doublereal(A(1,1));
	}
	else
	{
		doublereal* pB = B.Pointer();
		const double* pA = A.Pointer();
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB++ = doublereal(*pA++);
		*pB   = doublereal(*pA);
	}	
}

void ABAQUS_BaseT::ABAQUS_to_dSymMatrixT(const doublereal* pA,
	dSymMatrixT& B, bool convert_shear) const
{
	double shear_factor = (convert_shear) ? 0.5 : 1.0;
	double* pB = B.Pointer();
	if (B.Rows() == 2)
	{
		*pB++ = double(pA[0]); // 11
		*pB++ = double(pA[1]); // 22
		*pB   = double(pA[3])*shear_factor; // 12
	}
	else
	{
		*pB++ = double(pA[0]); // 11
		*pB++ = double(pA[1]); // 22
		*pB++ = double(pA[2]); // 33
		*pB++ = double(pA[5])*shear_factor; // 23
		*pB++ = double(pA[4])*shear_factor; // 13
		*pB   = double(pA[3])*shear_factor; // 12
	}
}

void ABAQUS_BaseT::dSymMatrixT_to_ABAQUS(const dSymMatrixT& A,
	doublereal* pB, bool convert_shear) const
{
	doublereal shear_factor = (convert_shear) ? 2.0 : 1.0;
	const double* pA = A.Pointer();
	if (A.Rows() == 2)
	{
		*pB++ = doublereal(pA[0]); // 11
		*pB++ = doublereal(pA[1]); // 22
		*pB++;                     // 33 - leave unchanged
		*pB   = doublereal(pA[2])*shear_factor; // 12
	}
	else
	{
		*pB++ = doublereal(pA[0]); // 11
		*pB++ = doublereal(pA[1]); // 22
		*pB++ = doublereal(pA[2]); // 33
		*pB++ = doublereal(pA[5])*shear_factor; // 12
		*pB++ = doublereal(pA[4])*shear_factor; // 13
		*pB   = doublereal(pA[3])*shear_factor; // 23
	}
}

/***********************************************************************
* Private
***********************************************************************/

/* read ABAQUS format input - reads until first line
* 	not containing a comment of keyword
*
* looks for kewords:
*    *MATERIAL
*    *USER MATERIAL
*    *DEPVAR
*/
void ABAQUS_BaseT::Read_ABAQUS_Input(ifstreamT& in, StringT& name, nArrayT<doublereal>& properties,
	double& density, integer& nstatv, bool& nonsym) const
{
	const char caller[] = "ABAQUS_BaseT::Read_ABAQUS_Input";

	/* disable comment skipping */
	int skip = in.skip_comments();
	char marker;
	if (skip)
	{
		marker = in.comment_marker();
		in.clear_marker();
	}

	/* required items */
	bool found_MATERIAL = false;
	bool found_USER     = false;

// when to stop
//    need at least material name
//    allow for no depvar and no constants
//    allow *DENSITY, *EXPANSION, *DAMPING <- mass damping

	/* advance to first keyword */
	bool found_keyword = Next_ABAQUS_Keyword(in);
	while (found_keyword)
	{
		StringT next_word;
		Read_ABAQUS_Word(in, next_word);
	
		/* first keyword */
		if (!found_MATERIAL)
		{
			if (next_word == "MATERIAL")
			{
				found_MATERIAL = true;
				/* scan for parameters */
				while (Skip_ABAQUS_Symbol(in, ','))
				{
					StringT next_parameter;
					Read_ABAQUS_Word(in, next_parameter);
					if (next_parameter == "NAME")
					{
						/* skip '=' */
						if (!Skip_ABAQUS_Symbol(in, '=')) ExceptionT::BadInputValue(caller);
						
						/* read name */
						Read_ABAQUS_Word(in, name, false);
					}
					else
						cout << "\n " << caller << ": skipping parameter:" << next_parameter << endl;
				}
			}
			else
				ExceptionT::BadInputValue(caller, "first keyword must be *MATERIAL not *%s",
					next_word.Pointer());
		}
		/* other keywords */
		else if (next_word == "USER")
		{
			nonsym = false;
			int nprops = 0;
			Read_ABAQUS_Word(in, next_word);
			if (next_word == "MATERIAL")
			{
				found_USER = true;
				/* scan for parameters */
				while (Skip_ABAQUS_Symbol(in, ','))
				{
					StringT next_parameter;
					Read_ABAQUS_Word(in, next_parameter);
					if (next_parameter == "CONSTANTS")
					{
						/* skip '=' */
						if (!Skip_ABAQUS_Symbol(in, '=')) ExceptionT::BadInputValue(caller);

						nprops = -1;
						in >> nprops;
						if (nprops < 0)
							ExceptionT::BadInputValue(caller, "error reading %s: %d",
								next_parameter.Pointer(), nprops);
					}
					else if (next_parameter == "UNSYMM")
						nonsym = true;
					else
						cout << "\n ABAQUS_BaseT::Read_ABAQUS_Input: skipping parameter:"
						     << next_parameter << endl;
				}
				
				/* read properties */
				properties.Dimension(nprops);
				in.clear_line();
				Skip_ABAQUS_Comments(in);
				for (int i = 0; i < nprops && in.good(); i++)
				{	
					in >> properties[i];
					Skip_ABAQUS_Symbol(in, ',');
					Skip_ABAQUS_Comments(in);
				}

				/* stream error */
				if (!in.good())
					ExceptionT::BadInputValue(caller, "error readind properties");
			}			
			else
				ExceptionT::BadInputValue(caller, "expecting MATERIAL after keyword *USER not %s",
					next_word.Pointer());
		}
		else if (next_word == "DEPVAR")
		{
			nstatv = -1;
			in >> nstatv;
			if (nstatv < 0)
				ExceptionT::BadInputValue(caller, "error reading %s", next_word.Pointer());

			/* clear trailing comma */
			Skip_ABAQUS_Symbol(in, ',');
		}
		else if (next_word == "DENSITY")
		{
			density = -1.0;
			in >> density;
			if (density < 0.0)
				ExceptionT::BadInputValue(caller, caller, "bad %s: %g", next_word.Pointer(), density);

			/* clear trailing comma */
			Skip_ABAQUS_Symbol(in, ',');
		}
		/* skipping all others */
		else
			cout << "\n ABAQUS_BaseT::Read_ABAQUS_Input: skipping keyword *"
			     << next_word << endl;

		/* next keyword */
		found_keyword = Next_ABAQUS_Keyword(in);
	}

	/* check all data found */
	if (!found_USER || !found_MATERIAL)
	{
		if (!found_MATERIAL)
			ExceptionT::BadInputValue(caller, "missing keyword: *MATERIAL");

		if (!found_USER)
			ExceptionT::BadInputValue(caller, "missing keyword: *USER MATERIAL");
	}

	/* restore comment skipping */
	if (skip) in.set_marker(marker);
}

/* advance to next non-comment line */
bool ABAQUS_BaseT::Next_ABAQUS_Keyword(ifstreamT& in) const
{
//rules: (1) keyword lines have a '*' in column 1
//       (2) comment lines have a '*' in column 1 and 2

	while (in.good())
	{
		/* advance to next '*' or alpanumeric character */
		char c = in.next_char();
		while (c != '*' && !isgraph(c) && in.good())
		{
			in.get(c);	
			c = in.next_char();
		}

		/* stream error */
		if (!in.good()) break;

		/* comment or keyword? */
		if (c == '*')
		{
			/* next */
			in.get(c);
			c = in.next_char();
			
			/* comment */
			if (c == '*')
				in.clear_line();
			/* keyword */
			else
				return true;
		}
		/* other alphanumeric character */
		else
			return false;
	}
	
	/* stream error */
	in.clear();
	return false;
}

void ABAQUS_BaseT::Skip_ABAQUS_Comments(ifstreamT& in) const
{
	while (in.good())
	{
		/* advance to next '*' or alpanumeric character */
		char c = in.next_char();
		while (c != '*' && !isgraph(c) && in.good())
		{
			in.get(c);	
			c = in.next_char();
		}

		/* stream error */
		if (!in.good()) break;
	
		/* comment or keyword? */
		if (c == '*')
		{
			in.get(c);
			c = in.next_char();
			
			/* comment */
			if (c == '*')
				in.clear_line();
			/* keyword */
			else {
				/* return the comment marker */
				in.putback('*');
				return;
			}
		}
		/* other alphanumeric character */
		else
			return;
	}

	/* stream error */
	in.clear();
}

bool ABAQUS_BaseT::Skip_ABAQUS_Symbol(ifstreamT& in, char c) const
{
	/* advance stream */
	char a = in.next_char();

	if (!in.good())
	{
		in.clear();
		return false;
	}
	else if (a == c)
	{
		in.get(a);
		return true;
	}
	else
		return false;
}

void ABAQUS_BaseT::Read_ABAQUS_Word(ifstreamT& in, StringT& word, bool to_upper) const
{
	const int max_word_length = 255;
	char tmp[max_word_length + 1];
	
	int count = 0;
	char c;
	/* skip whitespace to next word */
	c = in.next_char();
	
	in.get(c);
	while (count < max_word_length && (isalnum(c) || c == '-'))
	{
		tmp[count++] = c;
		in.get(c);
	}
	in.putback(c);
	tmp[count] = '\0';

	/* copy in */
	word = tmp;
	if (to_upper) word.ToUpper();
}

#endif /* __F2C__ */
