/* $Id: Anisotropic2DT.cpp,v 1.1.1.1 2001-01-29 08:20:25 paklein Exp $ */
/* created: paklein (06/11/1997)                                          */
/* Base class for 2D anisotropic materials                                */

#include "Anisotropic2DT.h"

#include <iostream.h>

#include "Rotate2DT.h"
#include "fstreamT.h"

/* constructors */
Anisotropic2DT::Anisotropic2DT(ifstreamT& in): fRotator(NULL)
{
	double theta12;
	
	in >> fIsRotated;
	if (fIsRotated != 0 && fIsRotated != 1) throw eBadInputValue;
	
	in >> theta12;	/* read in degrees */
		
	/* double check angle, but construct regardless */
	if (!fIsRotated) theta12 = 0.0;	
	
	fRotator = new Rotate2DT(theta12);		
	if (!fRotator) throw eOutOfMemory;
}

/* destructor */
Anisotropic2DT::~Anisotropic2DT(void)
{
	delete fRotator;
}

/* I/O functions */
void Anisotropic2DT::Print(ostream& out) const
{
	out << " Rotation with respect to global axes (rad). . . = ";
	out << fRotator->Angle() << '\n';
}

/*************************************************************************
* Protected
*************************************************************************/

/* transformation tensor */
const dMatrixT& Anisotropic2DT::Q(void) const
{
	return fRotator->Q();
}

/* return a reference to the transformed vector.  Note, returns a
* references to the argument if !fIsRotated */
const dArrayT& Anisotropic2DT::TransformIn(const dArrayT& vector)
{
	/* tranform into material coordinates */
	if (fIsRotated)
		return fRotator->RotateVectorIn(vector);
	else
		return vector;
}

const dArrayT& Anisotropic2DT::TransformOut(const dArrayT& vector)
{
	/* tranform out of material coordinates */
	if (fIsRotated)
		return fRotator->RotateVectorOut(vector);
	else
		return vector;
}	

/* return a reference to the transformed reduced index strain
* vector.  Note, returns a references to strain if !fRotator */
const dSymMatrixT& Anisotropic2DT::TransformIn(const dSymMatrixT& redmat)
{
	/* tranform into material coordinates */
	if (fIsRotated)
		return fRotator->RotateRedMatIn(redmat);
	else
		return redmat;
}

const dSymMatrixT& Anisotropic2DT::TransformOut(const dSymMatrixT& redmat)
{
	/* tranform out of material coordinates */
	if (fIsRotated)
		return fRotator->RotateRedMatOut(redmat);
	else
		return redmat;
}

/* 4th rank tensor tranformation - use for cases where moduli are constant
* and could therefore be stored in their transformed state */
void Anisotropic2DT::TransformIn(dMatrixT& redtensor)
{
	if (fIsRotated) fRotator->RotateRedTensorIn(redtensor);
}

void Anisotropic2DT::TransformOut(dMatrixT& redtensor)
{
	if (fIsRotated) fRotator->RotateRedTensorOut(redtensor);
}
