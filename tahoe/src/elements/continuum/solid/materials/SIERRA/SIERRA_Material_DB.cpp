/* $Id: SIERRA_Material_DB.cpp,v 1.5 2004-07-29 18:33:02 paklein Exp $ */
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
	for (int i = 0; i < fMaterialData.Size(); i++)
		delete tmp[i];	

	/* collect all pointers and free */
	ArrayT<C1FunctionT*> func;
	fFunctionEval.Ascending(func);
	for (int i = 0; i < fFunctionEval.Size(); i++)
		delete func[i];	
}
