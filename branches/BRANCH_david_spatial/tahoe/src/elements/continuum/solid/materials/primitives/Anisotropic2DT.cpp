/* $Id: Anisotropic2DT.cpp,v 1.6 2004-07-15 08:29:19 paklein Exp $ */
/* created: paklein (06/11/1997) */
#include "Anisotropic2DT.h"

#include <iostream.h>
#include "ifstreamT.h"
#include "Rotate2DT.h"

using namespace Tahoe;

/* constructors */
Anisotropic2DT::Anisotropic2DT(ifstreamT& in): 
	fRotator(NULL)
{
	int is_rotated;
	in >> is_rotated;
	if (is_rotated != 0 && is_rotated != 1) throw ExceptionT::kBadInputValue;
	
	double theta12;
	in >> theta12;	/* read in degrees */
		
	/* double check angle, but construct regardless */
	if (!is_rotated) theta12 = 0.0;	
	
	/* set rotator */
	fRotator = new Rotate2DT;
	SetRotation(theta12);
}

Anisotropic2DT::Anisotropic2DT(void):
	fRotator(NULL)
{
	/* set rotator */
	fRotator = new Rotate2DT;
}

/* set the rotation angle */
void Anisotropic2DT::SetRotation(double angle)
{
	/* construct new rotator */
	fRotator->SetAngle(angle);
}

/* destructor */
Anisotropic2DT::~Anisotropic2DT(void) { delete fRotator; }

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
	if (fabs(fRotator->Angle()) > kSmall)
		return fRotator->RotateVectorIn(vector);
	else
		return vector;
}

const dArrayT& Anisotropic2DT::TransformOut(const dArrayT& vector)
{
	/* tranform out of material coordinates */
	if (fabs(fRotator->Angle()) > kSmall)
		return fRotator->RotateVectorOut(vector);
	else
		return vector;
}	

/* return a reference to the transformed reduced index strain
* vector.  Note, returns a references to strain if !fRotator */
const dSymMatrixT& Anisotropic2DT::TransformIn(const dSymMatrixT& redmat)
{
	/* tranform into material coordinates */
	if (fabs(fRotator->Angle()) > kSmall)
		return fRotator->RotateRedMatIn(redmat);
	else
		return redmat;
}

const dSymMatrixT& Anisotropic2DT::TransformOut(const dSymMatrixT& redmat)
{
	/* tranform out of material coordinates */
	if (fabs(fRotator->Angle()) > kSmall)
		return fRotator->RotateRedMatOut(redmat);
	else
		return redmat;
}

/* 4th rank tensor tranformation - use for cases where moduli are constant
* and could therefore be stored in their transformed state */
void Anisotropic2DT::TransformIn(dMatrixT& redtensor)
{
	if (fabs(fRotator->Angle()) > kSmall)
		fRotator->RotateRedTensorIn(redtensor);
}

void Anisotropic2DT::TransformOut(dMatrixT& redtensor)
{
	if (fabs(fRotator->Angle()) > kSmall)
		fRotator->RotateRedTensorOut(redtensor);
}
