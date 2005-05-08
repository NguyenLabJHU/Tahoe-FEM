/* $Id: SIERRA_Material_DB.h,v 1.6 2004-08-08 02:02:57 paklein Exp $ */
#ifndef _SIERRA_MAT_DB_H_
#define _SIERRA_MAT_DB_H_

/* direct members */
#include "MapT.h"
#include "StringT.h"
#include "C1FunctionT.h"

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

	/** return a pointer to the material data associated with the given ID */
	static SIERRA_Material_Data* Material(int id);
	/*@}*/

	/** \name analytical functions */
	/*@{*/
	/** initialize new function with the given parameters */
	static void InitFunction(const ParameterListT& params);

	/** evaluate the given function */
	static double Evaluate(const StringT& name, double arg);
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

	/** map of material ID to data card */
	MapT<int, SIERRA_Material_Data*> fMaterialDataByID;

	/** map of name to function */
	MapT<StringT, C1FunctionT*> fFunctionEval;

	/** singleton to store parameters for Sierra materials */
	static SIERRA_Material_DB* the_SIERRA_Material_DB;
};

/* inlines */
/* return a reference to the singleton */
inline SIERRA_Material_DB& SIERRA_Material_DB::the_DB(void)
{
	if (!the_SIERRA_Material_DB)
		ExceptionT::GeneralFail("SIERRA_Material_DB::the_DB", "DB not instantiated");
	return *the_SIERRA_Material_DB;
}

inline SIERRA_Material_Data* SIERRA_Material_DB::Material(const StringT& name)
{
	/* get from DB */
	return the_DB().fMaterialData[name];
}

inline SIERRA_Material_Data* SIERRA_Material_DB::Material(int id)
{
	/* get from DB */
	return the_DB().fMaterialDataByID[id];
}

} /* namespace Tahoe */

#endif /* _SIERRA_MAT_DB_H_ */
