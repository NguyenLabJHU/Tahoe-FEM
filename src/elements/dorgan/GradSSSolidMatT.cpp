/* $Id: GradSSSolidMatT.cpp,v 1.4 2003-10-08 21:04:46 rdorgan Exp $ */ 
#include "GradSSSolidMatT.h"
#include <iostream.h>
#include "GradSSMatSupportT.h"
#include "dSymMatrixT.h"

using namespace Tahoe;

/* constructor */
GradSSSolidMatT::GradSSSolidMatT(ifstreamT& in, const GradSSMatSupportT& support):
        SSSolidMatT(in, support),
        fGradSSMatSupport(support),

        fNumDOF_R(support.NumDOF_R()),
        fNumDOF_Total(support.NumDOF() + fNumDOF_R),

        fNumIP_R(support.NumIP_R())
{

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
double& GradSSSolidMatT::R(void)
{
        return fGradSSMatSupport.LinearR(); 
}

double& GradSSSolidMatT::R(int ip)
{
        return fGradSSMatSupport.LinearR(ip); 
}

/* isotropic hardening from the end of the previous time step */
double& GradSSSolidMatT::R_last(void)
{
        return fGradSSMatSupport.LinearR_last(); 
}

double& GradSSSolidMatT::R_last(int ip)
{
        return fGradSSMatSupport.LinearR_last(ip); 
}

/* Laplacian isotropic hardening */
double& GradSSSolidMatT::LaplacianR(void)
{
        return fGradSSMatSupport.LinearLaplacianR(); 
}

double& GradSSSolidMatT::LaplacianR(int ip)
{
        return fGradSSMatSupport.LinearLaplacianR(ip); 
}

/* Laplacian isotropic hardening from the end of the previous time step */
double& GradSSSolidMatT::LaplacianR_last(void)
{
        return fGradSSMatSupport.LinearLaplacianR_last(); 
}

double& GradSSSolidMatT::LaplacianR_last(int ip)
{
        return fGradSSMatSupport.LinearLaplacianR_last(ip); 
}

/* apply pre-conditions at the current time step */
void GradSSSolidMatT::InitStep(void)
{
        /* inherited */
        SSSolidMatT::InitStep();
}
