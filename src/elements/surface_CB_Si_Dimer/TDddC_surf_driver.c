#include <stdlib.h>
void TDddC_surf_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi);
void TDddC_surf_1(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_2(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_3(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_4(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_5(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_6(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_7(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_8(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_9(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_10(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_11(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void TDddC_surf_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi) {
	double* z = (double*) malloc(21123*sizeof(double));
	TDddC_surf_1(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_2(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_3(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_4(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_5(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_6(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_7(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_8(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_9(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_10(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	TDddC_surf_11(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	free(z);
}
