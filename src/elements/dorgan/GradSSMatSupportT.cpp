/* $Id: GradSSMatSupportT.cpp,v 1.6 2004-04-01 22:46:54 rdorgan Exp $ */ 
#include "GradSSMatSupportT.h"
#include "ElementsConfig.h"

using namespace Tahoe;

/* constructor */
GradSSMatSupportT::GradSSMatSupportT(int nsd, int ndof_disp, int ndof_field, int nip_disp, int nip_field):
	SSMatSupportT(nsd, ndof_disp, nip_disp),

	fField_List(NULL),
	fField_last_List(NULL),
	fLapField_List(NULL),
	fLapField_last_List(NULL),

	fNumDOF_Field(ndof_field),
	fNumDOF_Total(ndof_disp + ndof_field),

	fNumIP_Field(nip_field)
{

}
 
/* destructor */
GradSSMatSupportT::~GradSSMatSupportT(void)
{

}

/* set source for the isotropic hadening */
void GradSSMatSupportT::SetLinearField(const dArrayT* field_List)
{
	/* keep pointer */
	fField_List = field_List;
}

/** set source for the isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearField_last(const dArrayT* field_last_List)
{
	/* keep pointer */
	fField_last_List = field_last_List;
}

/* set source for the laplacian of isotropic hadening */
void GradSSMatSupportT::SetLinearLaplacianField(const dArrayT* lapfield_List)
{
	/* keep pointer */
	fLapField_List = lapfield_List;
}

/** set source for the laplacian of isotropic hardening from the end of the previous time step */
void GradSSMatSupportT::SetLinearLaplacianField_last(const dArrayT* lapfield_last_List)
{
	/* keep pointer */
	fLapField_last_List = lapfield_last_List;
}

/* set the element group pointer */
void GradSSMatSupportT::SetContinuumElement(const ContinuumElementT* p)
{
	/* inherited */
	SSMatSupportT::SetContinuumElement(p);
}
