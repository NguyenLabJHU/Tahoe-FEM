/* $Id: SIERRA_Material_DB.cpp,v 1.1 2003-03-05 02:27:52 paklein Exp $ */
#include "SIERRA_Material_DB.h"
#include "SIERRA_Material_Data.h"

/* constructor */
SIERRA_Material_DB::SIERRA_Material_DB(void)
{

}

/* destructor */
SIERRA_Material_DB::~SIERRA_Material_DB(void)
{
	/* collect all pointers and free */
	ArrayT<SIERRA_Material_Data*> tmp;
	fMaterialData.Ascending(tmp);
	for (int i = 0; i < fMaterialData.Size(); i++)
		delete tmp[i];		
}
