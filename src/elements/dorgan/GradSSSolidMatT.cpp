/* $Id: GradSSSolidMatT.cpp,v 1.6 2004-01-14 19:33:16 rdorgan Exp $ */ 
#include "GradSSSolidMatT.h"
#include <iostream.h>
#include "GradSSMatSupportT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatT::GradSSSolidMatT(ifstreamT& in, const GradSSMatSupportT& support):
        SSSolidMatT(in, support),
        fGradSSMatSupport(&support),

        fNumDOF_R(support.NumDOF_R()),
        fNumDOF_Total(support.NumDOF() + fNumDOF_R),

        fNumIP_R(support.NumIP_R())
{
        SetName("gradient_small_strain_solid_material");
}

GradSSSolidMatT::GradSSSolidMatT(void):
        fGradSSMatSupport(NULL)
{
        SetName("gradient_small_strain_solid_material");
}

/* destructor */
GradSSSolidMatT::~GradSSSolidMatT(void) { }

/* initialization */
void GradSSSolidMatT::Initialize(void)
{
        /* inherited */
        SSSolidMatT::Initialize();
}

/* I/O */
void GradSSSolidMatT::PrintName(ostream& out) const
{
        /* inherited */
        SSSolidMatT::PrintName(out);
	
        out << "    Gradient Enhanced\n";
}

/* isotropic hardening */
const double& GradSSSolidMatT::R(void) const
{
        return fGradSSMatSupport->LinearR(); 
}

const double& GradSSSolidMatT::R(int ip) const
{
        return fGradSSMatSupport->LinearR(ip); 
}

/* isotropic hardening from the end of the previous time step */
const double& GradSSSolidMatT::R_last(void) const
{
        return fGradSSMatSupport->LinearR_last(); 
}

const double& GradSSSolidMatT::R_last(int ip) const
{
        return fGradSSMatSupport->LinearR_last(ip); 
}

/* Laplacian isotropic hardening */
const double& GradSSSolidMatT::LaplacianR(void) const
{
        return fGradSSMatSupport->LinearLaplacianR(); 
}

const double& GradSSSolidMatT::LaplacianR(int ip) const
{
        return fGradSSMatSupport->LinearLaplacianR(ip); 
}

/* Laplacian isotropic hardening from the end of the previous time step */
const double& GradSSSolidMatT::LaplacianR_last(void) const
{
        return fGradSSMatSupport->LinearLaplacianR_last(); 
}

const double& GradSSSolidMatT::LaplacianR_last(int ip) const
{
        return fGradSSMatSupport->LinearLaplacianR_last(ip); 
}

/* apply pre-conditions at the current time step */
void GradSSSolidMatT::InitStep(void)
{
        /* inherited */
        SSSolidMatT::InitStep();
}
