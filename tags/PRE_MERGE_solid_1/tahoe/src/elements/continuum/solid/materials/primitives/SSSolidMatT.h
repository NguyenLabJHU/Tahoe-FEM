/* $Id: SSSolidMatT.h,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/09/1997)                                          */
/* Defines the interface for elastic continuum materials.                 */

#ifndef _SS_STRUCT_MAT_T_H_
#define _SS_STRUCT_MAT_T_H_

/* base class */
#include "ContinuumT.h"
#include "StructuralMaterialT.h"

/* forward declarations */
class ShapeFunctionT;

class SSSolidMatT: protected ContinuumT, public StructuralMaterialT
{
public:

	/* constructor */
	SSSolidMatT(ifstreamT& in, const ElasticT& element);

	/* I/O functions */
	virtual void PrintName(ostream& out) const;

	/* required parameter flags */
	virtual bool NeedDisp(void) const;

	/* the shape functions */
	const ShapeFunctionT& ShapeFunction(void) const;
		
	/* elastic strain */
	const dSymMatrixT& e(void);

	/* material description */
	virtual const dMatrixT& C_IJKL(void);  // material tangent moduli
	virtual const dSymMatrixT& S_IJ(void); // PK2 stress
		// since infinitesimal strain materials don't make this
		// distinction, these functions just call the respective
		// c_ijkl and s_ij. this means 2 virtual function calls
		// per call. for efficiency, each derived SSSolidMatT should
		// overload these.

	/* apply pre-conditions at the current time step */
	virtual void InitStep(void);

protected:

	/* return the acoustical tensor and wave speeds */
	virtual const dSymMatrixT& AcousticalTensor(const dArrayT& normal);

private:

	/* set the internal thermal strain */
	virtual bool SetThermalStrain(dSymMatrixT& thermal_strain);

	/* acoustical tensor routines */
	void Q_2D(const dMatrixT& c_ijkl, const dArrayT& n, dSymMatrixT& Q) const;
	void Q_3D(const dMatrixT& c_ijkl, const dArrayT& n, dSymMatrixT& Q) const;

private:

	/* shape functions */
	const ShapeFunctionT& fShapes;
	
	/* nodal displacements */
	const LocalArrayT& fLocDisp;

	/* work space */
	dSymMatrixT	fStrainTemp; // elastic strain (w/o thermal)
	dSymMatrixT fQ;          // return value
	dMatrixT fGradU;         // displacement gradient matrix

	/* thermal strain: e_elastic = e_total - e_thermal */
	bool        fHasThermalStrain;
	dSymMatrixT fThermalStrain;

};

/* inlines */
inline const ShapeFunctionT& SSSolidMatT::ShapeFunction(void) const { return fShapes; }

#endif /* _SS_STRUCT_MAT_T_H_ */
