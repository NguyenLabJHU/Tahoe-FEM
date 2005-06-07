/* $Id: ParameterUtils.cpp,v 1.11 2005-06-07 07:29:53 paklein Exp $ */
#include "ParameterUtils.h"
#include <ctype.h>

using namespace Tahoe;

/**********************************************************************
 * IntegerListT implementation
 **********************************************************************/

IntegerListT::IntegerListT(const StringT& name):
	NamedListT<IntegerParameterT>(name)
{

}

IntegerListT::IntegerListT(void):
	NamedListT<IntegerParameterT>("IntegerList")
{

}

/**********************************************************************
 * DoubleListT implementation
 **********************************************************************/

/* constructors */
DoubleListT::DoubleListT(const StringT& name):
	NamedListT<DoubleParameterT>(name)
{

}

DoubleListT::DoubleListT(void):
	NamedListT<DoubleParameterT>("DoubleList")
{

}

/**********************************************************************
 * StringListT implementation
 **********************************************************************/

/* constructors */
StringListT::StringListT(const StringT& name):
	NamedListT<StringParameterT>(name)
{

}

StringListT::StringListT(void):
	NamedListT<StringParameterT>("StringList")
{

}

/* extract string parameters to an array */
void StringListT::Extract(const ParameterListT& list, ArrayT<StringT>& values)
{
	values.Dimension(list.NumLists("String"));
	for (int i = 0; i < values.Length(); i++)
		values[i] = list.GetList("String", i).GetParameter("value");
}

/**********************************************************************
 * IntegerT implementation
 **********************************************************************/

/* constructors */
IntegerParameterT::IntegerParameterT(void):
	NamedParameterT<ParameterT::Integer>("Integer")
{
	fValue = 0;
}

IntegerParameterT::IntegerParameterT(const StringT& name):
	NamedParameterT<ParameterT::Integer>(name)
{
	fValue = 0;
}

/**********************************************************************
 * DoubleT implementation
 **********************************************************************/

/* constructors */
DoubleParameterT::DoubleParameterT(void):
	NamedParameterT<ParameterT::Double>("Double")
{
	fValue = 0.0;
}

DoubleParameterT::DoubleParameterT(const StringT& name):
	NamedParameterT<ParameterT::Double>(name)
{
	fValue = 0.0;
}

/**********************************************************************
 * DoubleT implementation
 **********************************************************************/

/* constructors */
StringParameterT::StringParameterT(void):
	NamedParameterT<ParameterT::Word>("String")
{

}

StringParameterT::StringParameterT(const StringT& name):
	NamedParameterT<ParameterT::Word>(name)
{

}

/**********************************************************************
 * VectorParameterT implementation
 **********************************************************************/

VectorParameterT::VectorParameterT(const StringT& name, int dim, char variable):
	ParameterInterfaceT(name),
	fVariable(variable),
	fVector(dim)
{
	fVector = 0.0;
}

VectorParameterT::VectorParameterT(int dim, char variable):
	ParameterInterfaceT("vector"),
	fVariable(variable),
	fVector(dim)
{
	fVector = 0.0;
}

/* construct extracting length from the name */
VectorParameterT::VectorParameterT(const StringT& name_N, char variable):
	ParameterInterfaceT(name_N),
	fVariable(variable)
{
	const char caller[] = "VectorParameterT::VectorParameterT";
	const char msg[] = "could not extract length from \"%s\" in \"%s\"";

	/* resolve length */
	StringT suffix;
	suffix.Suffix(name_N, '_');
	if (suffix.StringLength() < 2 || !isdigit(suffix[1]))
		ExceptionT::GeneralFail(caller, msg, suffix, name_N.Pointer());
	int length = -1;
	length = atoi(suffix.Pointer(1));
	if (length < 0) ExceptionT::GeneralFail(caller, msg, suffix, name_N.Pointer());

	/* initialize */
	fVector.Dimension(length);
	fVector = 0.0;
}

/* describe the parameters needed by the interface */
void VectorParameterT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* define components */
	for (int i = 0; i < fVector.Length(); i++) {
		StringT v = "v_";
		v[0] = fVariable;
		v.Append(i+1);
		ParameterT v_i(ParameterT::Double, v);
		v_i.SetDefault(0.0);
		list.AddParameter(v_i);
	}
}

/* accept parameter list */
void VectorParameterT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	/* clear */
	fVector = 0.0;
	const ArrayT<ParameterT>& parameters = list.Parameters();
	for (int i = 0; i < parameters.Length(); i++) {
		const StringT& name = parameters[i].Name();
		if (name.StringLength() > 2 && name[0] == fVariable && name[1] == '_') {
			int component = atoi(name.Pointer(2)) - 1;
			if (component < 0 || component >= fVector.Length())
				ExceptionT::OutOfRange("VectorParameterT::TakeParameterList",
					"component \"%s\" is out of range {1,%d}", name.Pointer(), fVector.Length()+1);
			fVector[component] = parameters[i];
		}
	}
}

/* extract parameters to a dArrayT */
void VectorParameterT::Extract(const ParameterListT& list, dArrayT& array, char variable)
{
	VectorParameterT vec_param(list.Name(), variable);
	vec_param.TakeParameterList(list);
	array = vec_param;
}

/**********************************************************************
 * MatrixParameterT implementation
 **********************************************************************/

MatrixParameterT::MatrixParameterT(const StringT& name, int row, int col, char variable):
	ParameterInterfaceT(name),
	fVariable(variable),
	fMatrix(row, col),
	fCopySymmetric(false)
{
	fMatrix = 0.0;
}

MatrixParameterT::MatrixParameterT(int row, int col, char variable):
	ParameterInterfaceT("matrix"),
	fVariable(variable),
	fMatrix(row, col),
	fCopySymmetric(false)
{
	fMatrix = 0.0;
}

/* construct extracting dimensions from the name */
MatrixParameterT::MatrixParameterT(const StringT& name_NxM, char variable):
	ParameterInterfaceT(name_NxM),
	fVariable(variable),
	fCopySymmetric(false)
{
	const char caller[] = "MatrixParameterT::MatrixParameterT";
	const char msg[] = "could not extract %s dimensions from \"%s\" in \"%s\"";

	/* resolve suffix */
	StringT suffix;
	suffix.Suffix(name_NxM, '_');
	if (suffix.StringLength() < 4)
		ExceptionT::GeneralFail(caller, msg, "matrix", suffix, name_NxM.Pointer());
	
	/* resolve column dimensions */
	StringT num;
	num.Suffix(suffix, 'x');
	if (num.StringLength() < 2 || !isdigit(num[1]))
		ExceptionT::GeneralFail(caller, msg, "col", num.Pointer(), name_NxM.Pointer());
	int col = -1;
	col = atoi(num.Pointer(1));
	if (col < 0)
		ExceptionT::GeneralFail(caller, msg, "col", num.Pointer(), name_NxM.Pointer());
	
	/* resolve row dimensions */
	suffix.Root('x');
	if (suffix.StringLength() < 2 || !isdigit(suffix[1]))
		ExceptionT::GeneralFail(caller, msg, "row", suffix.Pointer(), name_NxM.Pointer());
	int row = -1;
	row = atoi(suffix.Pointer(1));
	if (row < 0)
		ExceptionT::GeneralFail(caller, msg, "row", suffix.Pointer(), name_NxM.Pointer());
	
	/* initialize */
	fMatrix.Dimension(row, col);
	fMatrix = 0.0;
}

/* describe the parameters needed by the interface */
void MatrixParameterT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* add flag */
	ParameterT copy_symmmetric(fCopySymmetric, "copy_symmetric");
	copy_symmmetric.SetDefault(fCopySymmetric);
	list.AddParameter(copy_symmmetric);

	/* define components */
	for (int i = 0; i < fMatrix.Cols(); i++)
		for (int j = 0; j < fMatrix.Rows(); j++) {
			StringT A = "A_";
			A[0] = fVariable;
			A.Append(j+1);
			A.Append("_", i+1);
			ParameterT A_ji(ParameterT::Double, A);
			A_ji.SetDefault(0.0);
			list.AddParameter(A_ji);
		}
}

/* accept parameter list */
void MatrixParameterT::TakeParameterList(const ParameterListT& list)
{
	/* inherited */
	ParameterInterfaceT::TakeParameterList(list);

	const char caller[] = "MatrixParameterT::TakeParameterList";
	const char msg[] = "%s of \"%s\" is out of range {1,%d}";

	/* clear */
	fMatrix = 0.0;
	StringT num, buffer;
	const ArrayT<ParameterT>& parameters = list.Parameters();
	for (int i = 0; i < parameters.Length(); i++) {
		const StringT& name = parameters[i].Name();
		if (name.StringLength() > 4 && name[0] == fVariable && name[1] == '_') {
			buffer = name;
			num.Suffix(buffer, '_');
			int col = atoi(num.Pointer(1)) - 1;
			buffer.Root('_');
			num.Suffix(buffer, '_');
			int row = atoi(num.Pointer(1)) - 1;
			
			/* checks */
			if (row < 0 && row >= fMatrix.Rows())
				ExceptionT::OutOfRange(caller, msg, "row", name.Pointer(), fMatrix.Rows()+1);
			if (col < 0 && col >= fMatrix.Cols())
				ExceptionT::OutOfRange(caller, msg, "col", name.Pointer(), fMatrix.Cols()+1);

			fMatrix(row, col) = parameters[i];
		}
	}

	/* copy symmetric */
	bool copy_symmetric = list.GetParameter("copy_symmetric");
	if (copy_symmetric && fMatrix.Rows() == fMatrix.Cols()) {
		for (int i = 0; i < fMatrix.Rows(); i++)
			for (int j = i+1; j < fMatrix.Rows(); j++) {
				double& a_ij = fMatrix(i,j);
				double& a_ji = fMatrix(j,i);
				if (fabs(a_ij) > kSmall && fabs(a_ji) > kSmall)
					ExceptionT::GeneralFail(caller, "copy symmetric with nonzero {%d,%d} and {%d,%d} entries",
						i+1, j+1, j+1, i+1);
				else if (fabs(a_ij) > kSmall)
					a_ji = a_ij;
				else /* a_ji > kSmall */
					a_ij = a_ji;				
			}
	}
}

/* extract parameters to a dArrayT */
void MatrixParameterT::Extract(const ParameterListT& list, dMatrixT& matrix, char variable)
{
	MatrixParameterT mat_param(list.Name(), variable);
	mat_param.TakeParameterList(list);
	matrix = mat_param;
}
