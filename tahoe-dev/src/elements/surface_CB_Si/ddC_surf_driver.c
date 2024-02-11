#include <stdlib.h>
void ddC_surf_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi);
void ddC_surf_1(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void ddC_surf_2(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void ddC_surf_3(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void ddC_surf_4(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void ddC_surf_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi) {
	double* z = (double*) malloc(8207*sizeof(double));
	ddC_surf_1(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	ddC_surf_2(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	ddC_surf_3(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	ddC_surf_4(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	free(z);
}
