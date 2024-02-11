#include <stdlib.h>
void ddC_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi);
void ddC_1(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void ddC_2(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi,double* z);
void ddC_driver(const double* params,const double* Xsi,const double* Xa,const double* Ya,const double* Za,const double* Cmat,double* dCdC,double* dCdXsi) {
	double* z = (double*) malloc(3770*sizeof(double));
	ddC_1(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	ddC_2(params,Xsi,Xa,Ya,Za,Cmat,dCdC,dCdXsi,z);
	free(z);
}
