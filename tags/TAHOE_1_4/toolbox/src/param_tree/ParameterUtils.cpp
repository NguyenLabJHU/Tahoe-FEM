/* $Id: ParameterUtils.cpp,v 1.8 2004-05-21 19:45:01 paklein Exp $ */
#include "ParameterUtils.h"

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

VectorParameterT::VectorParameterT(const StringT& name, char variable, int dim):
	ParameterInterfaceT(name),
	fVariable(variable),
	fVector(dim)
{
	fVector = 0.0;
}

VectorParameterT::VectorParameterT(char variable, int dim):
	ParameterInterfaceT("vector"),
	fVariable(variable),
	fVector(dim)
{
	fVector = 0.0;
}

/* construct extracting length from the name */
VectorParameterT::VectorParameterT(const StringT& name_N):
	ParameterInterfaceT(name_N),
	fVariable('v')
{
	const char caller[] = "VectorParameterT::VectorParameterT";
	const char msg[] = "could not extract length from \"%s\"";

	/* resolve length */
	StringT suffix;
	suffix.Suffix(name_N, '_');
	if (suffix.StringLength() < 2) ExceptionT::GeneralFail(caller, msg, name_N.Pointer());
	int length = -1;
	length = atoi(suffix.Pointer(1));
	if (length < 0) ExceptionT::GeneralFail(caller, msg, name_N.Pointer());

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
		ParameterT v_i = ParameterT(ParameterT::Double, v);
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
void VectorParameterT::Extract(const ParameterListT& list, dArrayT& array)
{
	VectorParameterT vec_param(list.Name());
	vec_param.TakeParameterList(list);
	array = vec_param;
}

/**********************************************************************
 * MatrixParameterT implementation
 **********************************************************************/

MatrixParameterT::MatrixParameterT(const StringT& name, char variable, int row, int col):
	ParameterInterfaceT(name),
	fVariable(variable),
	fMatrix(row, col)
{
	fMatrix = 0.0;
}

MatrixParameterT::MatrixParameterT(char variable, int row, int col):
	ParameterInterfaceT("matrix"),
	fVariable(variable),
	fMatrix(row, col)
{
	fMatrix = 0.0;
}

/* construct extracting dimensions from the name */
MatrixParameterT::MatrixParameterT(const StringT& name_NxM):
	ParameterInterfaceT(name_NxM),
	fVariable('A')	
{
	const char caller[] = "MatrixParameterT::MatrixParameterT";
	const char msg[] = "could not extract %s dimensions from \"%s\"";

	/* resolve suffix */
	StringT suffix;
	suffix.Suffix(name_NxM, '_');
	if (suffix.StringLength() < 4)
		ExceptionT::GeneralFail(caller, msg, "matrix", name_NxM.Pointer());
	
	/* resolve column dimensions */
	StringT num;
	num.Suffix(suffix, 'x');
	if (num.StringLength() < 2)
		ExceptionT::GeneralFail(caller, msg, "col", num.Pointer());
	int col = -1;
	col = atoi(num.Pointer(1));
	if (col < 0) ExceptionT::GeneralFail(caller, msg, "col", num.Pointer());
	
	/* resolve row dimensions */
	suffix.Root('x');
	if (suffix.StringLength() < 2)
		ExceptionT::GeneralFail(caller, msg, "row", suffix.Pointer());
	int row = -1;
	row = atoi(suffix.Pointer(1));
	if (row < 0) ExceptionT::GeneralFail(caller, msg, "row", suffix.Pointer());
	
	/* initialize */
	fMatrix.Dimension(row, col);
	fMatrix = 0.0;
}

/* describe the parameters needed by the interface */
void MatrixParameterT::DefineParameters(ParameterListT& list) const
{
	/* inherited */
	ParameterInterfaceT::DefineParameters(list);

	/* define components */
	for (int i = 0; i < fMatrix.Cols(); i++)
		for (int j = 0; j < fMatrix.Rows(); j++) {
			StringT A = "A_";
			A[0] = fVariable;
			A.Append(j+1);
			A.Append("_", i+1);
			ParameterT A_ji = ParameterT(ParameterT::Double, A);
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
}

/* extract parameters to a dArrayT */
void MatrixParameterT::Extract(const ParameterListT& list, dMatrixT& matrix)
{
	MatrixParameterT mat_param(list.Name());
	mat_param.TakeParameterList(list);
	matrix = mat_param;
}
