/* $Id: SIERRA_Material_DB.h,v 1.2 2003-03-06 17:23:31 paklein Exp $ */
#ifndef _SIERRA_MAT_DB_H_
#define _SIERRA_MAT_DB_H_

/* direct members */
#include "MapT.h"
#include "StringT.h"

namespace Tahoe {

/* forward declarations */
class SIERRA_Material_Data;

/** singleton database for Sierra materials parameters */
class SIERRA_Material_DB
{
public:

	/** \name managing the singleton */
	/*@{*/
	/** instantiate the singleton */
	static void Create(void);

	/** delete the singleton */
	static void Delete(void);
	/*@}*/

	/** \name registration functions 
	 * SIERRA_Material_DB::InitMaterial must be used to initialize
	 * the material before it can be accessed. */
	/*@{*/
	/** initialze new material */
	static void InitMaterial(const StringT& name, int XML_command_id, int modulus_flag);

	/** return a pointer to the material data associated with the given name */
	static SIERRA_Material_Data* Material(const StringT& name);
	/*@}*/

	/** \name indicies of real values */
	/*@{*/
	/** register a real constant */
	static void AddRealIndex(const StringT& name, int index);

	/** return the given real value */
	static int RealIndex(const StringT& name);
	/*@}*/

private:

	/** constructor */
	SIERRA_Material_DB(void);

	/** destructor */
	~SIERRA_Material_DB(void);

	/** return a pointer to the singleton */
	static SIERRA_Material_DB& the_DB(void);

private:

	/** array of material data cards. Maps material names to
	 * pointers to material data cards. */
	MapT<StringT, SIERRA_Material_Data*> fMaterialData;

	/** indicies of real values */
	MapT<StringT, int> fRealConstants;

	/** singleton to store parameters for Sierra materials */
	static SIERRA_Material_DB* the_SIERRA_Material_DB;
};

/* inlines */
inline SIERRA_Material_Data* SIERRA_Material_DB::Material(const StringT& name)
{
	/* get from DB */
	return the_DB().fMaterialData[name];
}

/* return a reference to the singleton */
inline SIERRA_Material_DB& SIERRA_Material_DB::the_DB(void)
{
	if (!the_SIERRA_Material_DB)
		ExceptionT::GeneralFail("SIERRA_Material_DB::the_DB", "DB not instantiated");
	return *the_SIERRA_Material_DB;
}

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_DB_H_ */
