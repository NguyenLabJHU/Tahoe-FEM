/* $Id: SIERRA_Material_DB.cpp,v 1.2 2003-03-06 17:23:31 paklein Exp $ */
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
}

/* register a real constant */
void SIERRA_Material_DB::AddRealIndex(const StringT& name, int index)
{
	/* store in the DB */
	the_DB().fRealConstants.Insert(name, index);
}

/* return the given real value */
int SIERRA_Material_DB::RealIndex(const StringT& name)
{
	return (the_DB().fRealConstants)[name];
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
}
