/* $Id: SIERRA_Material_DB.cpp,v 1.6 2004-08-08 02:02:57 paklein Exp $ */
#include "SIERRA_Material_DB.h"
#include "SIERRA_Material_Data.h"

using namespace Tahoe;

/* static data */
SIERRA_Material_DB* SIERRA_Material_DB::the_SIERRA_Material_DB = NULL;

/* instantiate the singleton */
void SIERRA_Material_DB::Create(void)
{
	if (!the_SIERRA_Material_DB) 
		the_SIERRA_Material_DB = new SIERRA_Material_DB;
}

/* delete the singleton */
void SIERRA_Material_DB::Delete(void)
{
	delete the_SIERRA_Material_DB;
	the_SIERRA_Material_DB = NULL;
}

/* initialze new material */
void SIERRA_Material_DB::InitMaterial(const StringT& name, int XML_command_id, int modulus_flag)
{
	/* create material data card */
	SIERRA_Material_Data* data = new SIERRA_Material_Data(name, XML_command_id, modulus_flag);

	/* store in the DB */
	the_DB().fMaterialData.Insert(name, data);
	
	/* store by ID */
	the_DB().fMaterialDataByID.Insert(data->ID(), data);
}

/* initialize new function */
void SIERRA_Material_DB::InitFunction(const ParameterListT& params)
{
	const char caller[] = "SIERRA_Material_DB::InitFunction";
	
	/* check */
	if (params.Description() != "FUNCTION")
		ExceptionT::GeneralFail(caller, "expecting description \"FUNCTION\" not \"%s\"",
			params.Description().Pointer());
			
	/* function type */
	const StringT& function_type = params.GetParameter("TYPE");
	ParameterListT c1_params;
	if (function_type == "PIECEWISE_LINEAR")
	{
		/* translate parameters */
		c1_params.SetName("piecewise_linear");
		const ParameterListT& points = params.GetList("VALUES");
		int npts = points.NumLists("OrderedPair");
		for (int i = 0; i < npts; i++) /* copy order pairs */
			c1_params.AddList(points.GetList("OrderedPair", i));
	}
	else
		ExceptionT::GeneralFail(caller, "unrecognized function type \"%s\"",
			function_type.Pointer());

	/* construct function */
	C1FunctionT* c1 = C1FunctionT::New(c1_params.Name());
	if (!c1) 
		ExceptionT::GeneralFail(caller, "could not construct \"%s\"",
			c1_params.Name().Pointer());
	c1->TakeParameterList(c1_params);

	/* store in the DB */
	the_DB().fFunctionEval.Insert(params.Name(), c1);
}

/* evaluate the given function */
double SIERRA_Material_DB::Evaluate(const StringT& name, double arg)
{
	C1FunctionT* func = the_DB().fFunctionEval[name];
	return func->Function(arg);
}

/***********************************************************************
 * Private
 ***********************************************************************/

/* constructor */
SIERRA_Material_DB::SIERRA_Material_DB(void) { }

/* destructor */
SIERRA_Material_DB::~SIERRA_Material_DB(void)
{
	/* collect all pointers and free */
	ArrayT<SIERRA_Material_Data*> tmp;
	fMaterialData.Ascending(tmp);
	for (int i = 0; i < tmp.Length(); i++)
		delete tmp[i];	

	/* collect all pointers and free */
	ArrayT<C1FunctionT*> func;
	fFunctionEval.Ascending(func);
	for (int i = 0; i < func.Length(); i++)
		delete func[i];	
}
